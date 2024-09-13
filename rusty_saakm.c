#include "asm-generic/errno-base.h"
#include "linux/cpumask.h"
#include "linux/gfp_types.h"
#include "linux/kernel_stat.h"
#include "linux/percpu-defs.h"
#include "linux/printk.h"
#include "linux/sched.h"
#include "linux/sched/smt.h"
#include "linux/smp.h"
#include "linux/topology.h"
#include <linux/module.h>
#include <linux/slab.h>
#include <linux/proc_fs.h>
#include <linux/debugfs.h>
#include <linux/ipanema.h>
#include <linux/kthread.h>
#include <linux/delay.h>
#include <../kernel/sched/sched.h>

MODULE_AUTHOR("Raïssa Bouri, lip6");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("rusty saakmed");

#define SLICE_NS  20 * 1000000
#define MAX_DOMS 15
#define PICK_IDLE_CORE 1LLU<<0

/* for the tuner */
static struct task_struct *tuner_thread;

static int nr_doms;
static bool idle_enabled = true; //<=> scx_builtin_idle_enabled

/* ipanema stuff */
struct ipanema_policy *rusty_policy;
static const char *policy_name = "rusty\0";

/* representation of a core */
struct rusty_core {
    enum ipanema_core_state state;
    int id;
    struct ipanema_rq rq;

    /* rusty pcpu_ctx fields */
    // domains are round robined
    //this field will help go through them
    int dom_rr_cur; // <=> ptr sur next
    u32 dom_id;

    //stats for the tuner
    u64 user   ;  
    u64 nice   ;  
    u64 system ;  
    u64 idle   ;  
    u64 iowait ;  
    u64 irq    ;  
    u64 softirq;  
    u64 steal  ;  
};

DEFINE_PER_CPU(struct rusty_core, core);

/* representation of a process */
struct rusty_process {
    enum ipanema_state state;
    struct ipanema_rq *rq;
    struct task_struct *task;

    /* rusty task_ctx fields */
    cpumask_t * cpumask;
    cpumask_t * tmp_cpumask;
    struct topology_level *l; //<=> rusty's dom_ctx
    u64 dom_mask;//updated in task_pick_domain
    int dom_id;
    int slice;
    bool dispatch_local;
    bool all_cpus;
};

/* idle masks */
struct {
    cpumask_t cpu;
    cpumask_t smt;
} idle_masks;

/* tune inputs */
static cpumask_t * direct_greedy_cpumask;
/* kick_greedy is used in main.bpf.c's rusty_enqueue only in case the task's dispatch_local is false
 it only happens if rusty_select_cpu doesn't go through the direct label
 this doesn't happen in my code
*/
static cpumask_t * kick_greedy_cpumask;
static cpumask_t * all_cpumask;

/* dom_ctx 
in rusty this structure contains a struct bucket for the load-balancing

also, it's to round-robin the domains
*/
struct dom_ctx {
    u32 id;
    struct topology_level * topology;
};
struct dom_ctx *dom_data[MAX_DOMS];

/* rusty default parameters but the user can change them
they're in main.rs
*/
static bool direct_greedy_numa = false;
static int direct_greedy_under = 90 ;
static int kick_greedy_under = 100;

/* tuning part, see tuner.rs and the big note i wrote before rusty_routine at the end of this file*/
static void tuner_step(void) {

    int cpu, dom_cpu, dom_util = 0;
    struct topology_level *l;
    struct rusty_core *c;

    //for stats
    cpumask_t *dom_already_done = kzalloc(sizeof(cpumask_size()), GFP_KERNEL);
    struct kernel_cpustat kcpustat;
    u64 *cpustat;

    cpumask_clear(dom_already_done);

    //current cpustat
    u64 new_user, new_nice, new_system, new_idle, new_iowait, new_irq, new_softirq, new_steal;
    new_user = new_nice = new_system = new_idle = new_iowait = new_irq = new_softirq = new_steal = 0;

    //util cpustat
    bool enable_direct, enable_kick;
    enable_direct = enable_kick = false;

    u64 util_user, util_nice, util_system, util_idle, util_iowait, util_irq, util_softirq, util_steal;
    util_user = util_nice = util_system = util_idle = util_iowait = util_irq = util_softirq = util_steal = 0;

    u64 busy_use, total_use, util, dom_sum;
    busy_use = total_use = util = dom_sum = 0;

    for_each_possible_cpu(cpu) {
     
        //get domain mask
        l = per_cpu(topology_levels, cpu);
        
        //if cpu was already checked below then no need to do it again (i guess)
        if (cpumask_test_cpu(cpu, dom_already_done)) 
            continue;
        
        for_each_cpu(dom_cpu, &l->cores){ //TODO check if numa level ?

            //set cpu in dom_already_done 
            cpumask_set_cpu(dom_cpu, dom_already_done);

            //get new stats of CPUs in that domain
            cpustat = kcpustat.cpustat;
            kcpustat_cpu_fetch(&kcpustat, cpu);

            //to get prev stats stored in the rusty_proces struct
            c = &per_cpu(core, dom_cpu);

            //the new stats
            new_user    = div_u64(cpustat[CPUTIME_USER]   , NSEC_PER_SEC / USER_HZ);         
            new_nice    = div_u64(cpustat[CPUTIME_NICE]   , NSEC_PER_SEC / USER_HZ);
            new_system  = div_u64(cpustat[CPUTIME_SYSTEM] , NSEC_PER_SEC / USER_HZ);
            new_idle    = div_u64(cpustat[CPUTIME_IDLE]   , NSEC_PER_SEC / USER_HZ);
            new_iowait  = div_u64(cpustat[CPUTIME_IOWAIT] , NSEC_PER_SEC / USER_HZ);
            new_irq     = div_u64(cpustat[CPUTIME_IRQ]    , NSEC_PER_SEC / USER_HZ);
            new_softirq = div_u64(cpustat[CPUTIME_SOFTIRQ], NSEC_PER_SEC / USER_HZ);
            new_steal   = div_u64(cpustat[CPUTIME_STEAL]  , NSEC_PER_SEC / USER_HZ);

            //<=> sub_or_zero between cur and prev, that's a function in tuner.rs
            util_user = (new_user - c->user) < 0 ? 0 : (new_user - c->user) ;
            util_nice = (new_nice - c->nice) < 0 ? 0 : (new_nice - c->nice) ;
            util_system = (new_system - c->system) < 0 ? 0 : (new_system - c->system) ;
            util_idle = (new_idle - c->idle) < 0 ? 0 : (new_idle - c->idle) ;
            util_iowait = (new_iowait - c->iowait) < 0 ? 0 : (new_iowait - c->iowait) ;
            util_irq = (new_irq - c->irq) < 0 ? 0 : (new_irq - c->irq) ;
            util_softirq = (new_softirq - c->softirq) < 0 ? 0 : (new_softirq - c->softirq) ;
            util_steal = (new_steal - c->steal) < 0 ? 0 : (new_steal - c->steal) ;

            //update new stats in core struct
            c->user     = new_user ;
            c->nice     = new_nice ;
            c->system   = new_system ;
            c->idle     = new_idle;
            c->iowait   = new_iowait ;
            c->irq      = new_irq ;
            c->softirq  = new_softirq ;
            c->steal    = new_steal ;

            //tuner computations
            busy_use = util_user + util_system + util_nice + util_irq + util_softirq + util_steal;
            total_use = util_idle + busy_use + util_iowait ;

            if (total_use > 0) 
                util =  max((busy_use / total_use) , 1);
            else 
                util = 1;

            dom_util += util;

        }

        if (cpumask_weight(&l->cores) > 0 )
            util = dom_util / cpumask_weight(&l->cores);
        else 
            util = 0;
        
        //TODO check if floats okay
        enable_direct = direct_greedy_under > 0.99999 || util < direct_greedy_under;
        enable_kick = kick_greedy_under > 0.99999 || util < kick_greedy_under;
        
        if (enable_direct) 
            cpumask_or(direct_greedy_cpumask, direct_greedy_cpumask, &l->cores);

        if (enable_kick)
            cpumask_or(kick_greedy_cpumask , kick_greedy_cpumask, &l->cores);

    }
}

/* core func */

static enum ipanema_core_state sched_get_core_state (struct ipanema_policy *policy ,
                                                 struct core_event *e )
{
    return per_cpu(core, e->target).state;
}

static void sched_core_entry(struct ipanema_policy *policy, struct core_event *e)
{
 
    struct rusty_core *c = &per_cpu(core, e->target);
    c->state = IPANEMA_ACTIVE_CORE;
    cpumask_set_cpu(c->id, all_cpumask);

    if (sched_smt_active()) {
        cpumask_clear_cpu(c->id, &idle_masks.smt);
    }
    else {
        cpumask_clear_cpu(c->id, &idle_masks.cpu);
    }
}

static void sched_core_exit(struct ipanema_policy *policy, struct core_event *e)
{
    struct rusty_core *c = &per_cpu(core, e->target);
    c->state = IPANEMA_IDLE_CORE;
    cpumask_clear_cpu(c->idle, all_cpumask);

    if (sched_smt_active()) {
        cpumask_clear_cpu(c->id, &idle_masks.smt);
    }
    else {
        cpumask_clear_cpu(c->id, &idle_masks.cpu);
    }
}

static void sched_newly_idle(struct ipanema_policy *policy, struct core_event *e)
{
}

static void sched_enter_idle(struct ipanema_policy *policy, struct core_event *e)
{
    struct rusty_core *c = &per_cpu(core,e->target);
    per_cpu(core, e->target).state = IPANEMA_IDLE_CORE;
    if (sched_smt_active()) {
        cpumask_set_cpu(c->id, &idle_masks.smt);
    }
    else {
        cpumask_set_cpu(c->id, &idle_masks.cpu);
    }
}

static void sched_exit_idle(struct ipanema_policy *policy, struct core_event *e)
{
     struct rusty_core *c = &per_cpu(core,e->target);
    per_cpu(core, e->target).state = IPANEMA_ACTIVE_CORE;

    if (sched_smt_active()) {
        cpumask_clear_cpu(c->id, &idle_masks.smt);
    }
    else {
        cpumask_clear_cpu(c->id, &idle_masks.cpu);
    }
}

static void sched_balancing_select(struct ipanema_policy *policy, struct core_event *e)
{
}

/* ext.c's function */

static bool test_and_clear_cpu_idle(int cpu) 
{
    if (sched_smt_active()) {
        const struct cpumask *smt = cpu_smt_mask(cpu);

        if (cpumask_intersects(smt, &idle_masks.smt))
            cpumask_andnot(&idle_masks.smt, &idle_masks.smt, smt);

        else if (cpumask_test_cpu(cpu, &idle_masks.smt))
            cpumask_clear_cpu(cpu, &idle_masks.smt);
    }

    return cpumask_test_and_clear_cpu(cpu, &idle_masks.cpu);
}

static cpumask_t * get_idle_cpumask(void) {

    if (sched_smt_active()) 
        return &idle_masks.smt;
    else
        return &idle_masks.cpu;
}

static int pick_any_cpu(const struct cpumask *cpus_allowed, u64 flags)
{
    int cpu;
    /* didn't have time to look into the scx_builtin_idle_enabled but 
    *  it's likely so i copied it
    */

    if (idle_enabled) {
        cpu = pick_any_cpu(cpus_allowed, flags);
        if (cpu >= 0)
            return cpu;
    }

    cpu = cpumask_any_distribute(cpus_allowed);
    if (cpu < nr_cpu_ids)
        return cpu;
    else
        return -EBUSY;


}
/* just like main.bpf.c's try_sync_wakeup fun */
static s32 try_sync_wakeup (struct task_struct *p, struct rusty_process *rp, s32 prev_cpu) {

    int cpu;
    struct rusty_core *c = kzalloc(sizeof(struct rusty_core), GFP_KERNEL);
    struct dom_ctx *domc = kzalloc(sizeof(struct dom_ctx), GFP_KERNEL);
    cpumask_t *d_cpumask, *idle_cpumask;
    d_cpumask = kzalloc(sizeof(cpumask_size()), GFP_KERNEL);
    idle_cpumask = kzalloc(sizeof(cpumask_size()), GFP_KERNEL);

    bool share_llc, has_idle;

    cpu = smp_processor_id();
    c = &per_cpu(core,cpu);
    domc = dom_data[c->dom_id];
    d_cpumask = &domc->topology->cores;

    idle_cpumask = get_idle_cpumask();

    share_llc = cpumask_test_cpu(prev_cpu, (const struct cpumask *)d_cpumask);
    if (share_llc && test_and_clear_cpu_idle(prev_cpu)) {
        cpu = prev_cpu;
        goto err_out;
    }
   
    has_idle = cpumask_intersects((const struct cpumask *)d_cpumask, (const struct cpumask *)idle_cpumask);

    if (has_idle && cpumask_test_cpu(cpu, p->cpus_ptr) 
        && !(current->flags & PF_EXITING) && (rp->dom_id < MAX_DOMS)
        && ( c->rq.nr_tasks ==0)) {
           
            goto err_out;
    } 

    cpu = -ENOENT;

err_out:
    return cpu;
}

static s32 pick_idle_cpu (const struct cpumask *cpus_allowed, u64 flags)
{
    int cpu;

retry :
    if (sched_smt_active()) {
        cpu = cpumask_any_and_distribute(&idle_masks.smt, cpus_allowed);
    
        if (cpu < nr_cpu_ids)
            goto found;

        if (flags & PICK_IDLE_CORE)
            return -EBUSY;
    }
  
    cpu = cpumask_any_and_distribute(&idle_masks.cpu, cpus_allowed);
    if (cpu >= nr_cpu_ids){
        return -EBUSY;  
}
found :

    if (test_and_clear_cpu_idle(cpu)) {
        
        return cpu;
        }
    else
     goto retry;;
}

/* I put the same comments as Rusty as a landmark  */
static int rusty_select_cpu(struct rusty_process *rp, struct task_struct *p, int prev_cpu, int wake_flags) 
{
    const struct cpumask *idle_smtmask = kzalloc(sizeof(cpumask_size()), GFP_KERNEL);
    struct cpumask *p_cpumask =kzalloc(sizeof(cpumask_size()), GFP_KERNEL);

    s32 cpu;
    bool has_idle_cores, prev_domestic;

    p_cpumask = rp->cpumask;

    if (sched_smt_active()) 
        idle_smtmask = &idle_masks.smt;
    else
        idle_smtmask = &idle_masks.cpu;

    //to refresh the tune params
    tuner_step();

    if (p->nr_cpus_allowed == 1) {
        cpu = prev_cpu;
        goto direct;
    }
    
    if (wake_flags & WF_SYNC) {
        cpu = try_sync_wakeup(p, rp, prev_cpu);
        if (cpu >= 0)
            goto direct;
    }

    has_idle_cores = !cpumask_empty(idle_smtmask);

    //did task get pulled to foreing domain ?
    prev_domestic = cpumask_test_cpu(prev_cpu, (const struct cpumask *)p_cpumask);

    //if whole physical core is idle we keep prev_cpu
    if (prev_domestic) {
        if (cpumask_test_cpu(prev_cpu, idle_smtmask)
            && test_and_clear_cpu_idle(prev_cpu) ) {
     
            cpu = prev_cpu;
            goto direct;
        } 
    } 

    /*prev cpu won't do
    looking for an idle cpu we can directly dispatch the task to
    first : look for best idle domestic cpu
    and then move onto foreign
    */

    //if there's a domestic idle core then dispatch directly
    if (has_idle_cores) {
     
        cpu = pick_idle_cpu((const struct cpumask *)p_cpumask, PICK_IDLE_CORE); //TODO
        if (cpu >= 0) {
            goto direct;
        }
    }

    //actually if prev cpu was domestic and idle even if core wasn't then it might be okay
    if (prev_domestic && test_and_clear_cpu_idle(prev_cpu)) {
        cpu = prev_cpu;
        goto direct;
    }
 
    //if any domestic idle cpu dispatch directly
    cpu = pick_idle_cpu((const struct cpumask*)p_cpumask, 0); //TODO
    if (cpu >= 0)
        goto direct;

    /* domestic domain fully booked
    if there are CPUs which are idle and under utilized (direct_greedy) then ignore domain boundaries
    while still respecting numa's and push the task there
    */

    //first, try to find an idle core

    if (rp->all_cpus && direct_greedy_cpumask 
        && (!cpumask_empty((const struct cpumask *)direct_greedy_cpumask))) {

        cpumask_t *tmp_direct_greedy, *node_mask ;
        tmp_direct_greedy= kzalloc(sizeof(cpumask_size()), GFP_KERNEL);//, *tmp_cpumask = NULL;
        node_mask=kzalloc(sizeof(cpumask_size()), GFP_KERNEL);

        tmp_direct_greedy = direct_greedy_cpumask ;

        /*by default only look for an idle core in the current numa node 
        when looking for direct_greedy cpus outside of the current domain
        this is set by direct_greedy_numa which default value is false in main.rs
        */

        if (!direct_greedy_numa) {

            //they need the numa node mask
            //i can get it if topology-level->flag & DOMAIN_NUMA
    
            cpumask_and(rp->tmp_cpumask, 
                        (const struct cpumask *)node_mask,
                        (const struct cpumask *)tmp_direct_greedy);
            tmp_direct_greedy = rp->tmp_cpumask;
        }


        /*---------- note on this part -----------*/
        /* depending on whether has_idle_cores is true 
        * the flag passed to pick_idle_cpu (ext.c function) changes
        * when PICK_IDLE_CORE is passed there's no possibility to remain stuck in a loop
        * because if no cpu is found it'll return EBUSY
        */

        //try to find an idle core in the previous and then any domain
        if (has_idle_cores) {
            /* not fully done because i didn't add a direct_greedy_cpumask for each domain
            * but in rusty they first look into it ( = idle core in the previous dom)
            * then they look in the global direct_greedy_cpumask (= any dom)
            */

            if (direct_greedy_cpumask) {
                cpu = pick_idle_cpu((const struct cpumask*)tmp_direct_greedy,
                                    PICK_IDLE_CORE);
                if (cpu >= 0) {
                    goto direct;
                }
            }
        }

        //no idle core...is there any ide cpu ?
        //skipping the first part aswell because they look into the dom's direct_greedy_cpumask
        if (direct_greedy_cpumask) {
          
            cpu = pick_idle_cpu((const struct cpumask *)tmp_direct_greedy, 0);

            if (cpu >= 0) {
                goto direct;
            }
        }        
    }

    /*----------end of the note----------*/

    /* queue on domestic domain's dispatch queue even tho it might 
    * be on a different domain, this might lead to stalls :/ 
    */
    if (prev_domestic) 
        cpu = prev_cpu;
    else 
        cpu = pick_any_cpu((const struct cpumask *)p_cpumask, 0);

    return cpu;

direct : 
    rp->dispatch_local = true;
    return cpu;

}

/* just like Rusty */
static void task_pick_domain (struct rusty_process *rp, struct task_struct *p)
{

    int cpu = smp_processor_id();
    int  i;
    int first_dom = -1;
    struct rusty_core * c = &per_cpu(core, cpu);
    u32 dom;

    rp->dom_mask = 0;
    dom = c->dom_rr_cur++;

    for (i = 0; i < nr_doms ; i++) {
        dom = (dom + 1) % nr_doms;
        if (cpumask_intersects(p->cpus_ptr, &dom_data[dom]->topology->cores)) {
           
            rp->dom_mask |= 1LLU << dom;
            
            cpumask_and(rp->cpumask, (const struct cpumask *)&dom_data[dom]->topology->cores,
                                     p->cpus_ptr);
            
            if (first_dom < 0) {
                first_dom = dom;
                rp->dom_id = dom;
            }
        }
    }
}

static int sched_new_prepare (struct ipanema_policy *policy, struct process_event *e )
{
    int cpu;

    /* init process struct*/
    struct rusty_process *rp;
    struct task_struct *p = e->target;
    rp = kzalloc(sizeof(struct rusty_process), GFP_ATOMIC);
    if (!rp)
        return -1;

    
    rp->task = p;
    rp->rq = NULL;
    rp->dispatch_local = false;
    rp->state = IPANEMA_NOT_QUEUED;
    p->ipanema.policy_metadata = rp;

    rp->cpumask = kzalloc(sizeof(cpumask_size()), GFP_KERNEL);
    rp->tmp_cpumask = kzalloc(sizeof(cpumask_size()), GFP_KERNEL);
    cpumask_clear(rp->cpumask);
    cpumask_clear(rp->tmp_cpumask); //TODO

    /* init task's cpumasks */
    /* 
        - done by task_pick_and_set_domain(task_ctx, task_struct p,
                                        p->cpus-ptr, true) 
        
            dom_id = task_pick_domain()
            - task_pick_domain(task_ctx, p, p->cpus_ptr)
                * gets smp_processor_id
                * dom = next domain in round robined table 
                    for each domain
                        dom = index of next domain in rr table 
                        if (dom in p->cpus_ptr)
                            set dom in dom_mask
    */

    task_pick_domain(rp, p);

    //set rusty_process->all_cpus to true if all the cpus in all_cpumask are in p->cpus_ptr
    rp->all_cpus = cpumask_subset((const struct cpumask *)all_cpumask, p->cpus_ptr);

    cpu = rusty_select_cpu(rp, p, task_cpu(p), e->flags);
    return cpu;
}

/* helper function to enqueue a task                                        
   kinda the rusty equivalent 
   in rusty_enqueue
   -> check if lb_data requested a migration (skipped)
   -> check if dispatch_local true 
   (here it should be, in rusty it might be false in bad cases -> can lead to stalls)
        -> if true 
            dispatch task on local dsq <=> core's ipanema_rq with a time slice of slice_ns
            which is equal to SCX_SLICE_DFL = 20 * 1000000,	(20 ms) defined in linux/sched/ext.h
            call scx_bpf_dispatch from ext.c
*/
static void rusty_enqueue_task(struct rusty_process *rp, struct ipanema_rq *next_rq, int state)
{
    rp->state = state;
    rp->rq = next_rq;

    change_state(rp->task, state, task_cpu(rp->task), next_rq);
}

/*
    called right after ipanema_new_prepare
    to enqueue the task on the chosen cpu
*/
static void sched_new_place (struct ipanema_policy *policy, struct process_event *e)
{
    struct rusty_process *rp = policy_metadata(e->target);
    int cpu_dst = task_cpu(rp->task);
    
    //set slice
    rp->slice = SLICE_NS;

    rusty_enqueue_task(rp, &per_cpu(core, cpu_dst).rq, IPANEMA_READY);

    /* in ext.c they call rusty_runnable from enqueue_task */
    /* it has to do with load balancing so i skipped that part
        but basically the dom_ctx struct contains an array bucket[LB_LOAD_BUCKET] 
        and they do some computatoins on this
    */
}

static void sched_new_end (struct ipanema_policy *policy , struct process_event *e)
{ 
}

static void sched_tick (struct ipanema_policy *policy, struct process_event *e)
{  
    /* task_tick_scx in ext.c
        -> calls update_curr_scx : what ipanema does + this
    */
    struct rusty_process *rp = policy_metadata(e->target);
    u64 delta_exec;

    /* todo : now should be rq_clock_task but i cant use that 
    it takes a rq and they're not exported
    */
    delta_exec = ktime_get_ns() - e->target->se.exec_start;
    rp->slice -= (rp->slice < delta_exec) ? rp->slice : delta_exec;
    
    /* -> second part of their tick function */
    if (!rp->slice)
        rusty_enqueue_task(rp, &per_cpu(core, task_cpu(rp->task)).rq, IPANEMA_READY_TICK);
}

/* in ext.c yield_task_scx check if the scheduler has the yield function
    rusty doesn't so it just sets the task's slice to 0
*/
static void sched_yield (struct ipanema_policy *policy, struct process_event *e)
{
    struct rusty_process *rp = policy_metadata(e->target);

    rusty_enqueue_task(rp, &per_cpu(core, task_cpu(rp->task)).rq, IPANEMA_READY);

    rp->slice = 0;
}

static void sched_block(struct ipanema_policy *p, struct process_event *e)
{
    struct rusty_process *rp = policy_metadata(e->target);

    rp->state = IPANEMA_BLOCKED;
    rp->rq = NULL;

    change_state(rp->task, IPANEMA_BLOCKED, task_cpu(rp->task), NULL);
}

static int sched_unblock_prepare(struct ipanema_policy *p, struct process_event *e)
{   
    return rusty_select_cpu(policy_metadata(e->target), e->target, task_cpu(e->target), 0);
}

static void sched_unblock_place(struct ipanema_policy *p, struct process_event *e)
{ 
    struct rusty_process *rp = policy_metadata(e->target);
    int cpu_dst = task_cpu(rp->task);

    rusty_enqueue_task(rp, &per_cpu(core, cpu_dst).rq, IPANEMA_READY);
}

static void sched_unblock_end(struct ipanema_policy *policy, struct process_event *e)
{
}

static void sched_terminate(struct ipanema_policy *policy, struct process_event *e)
{
   struct rusty_process *rp = policy_metadata(e->target);
    
    rp->state = IPANEMA_TERMINATED;
    rp->rq = NULL;

    change_state(rp->task, IPANEMA_TERMINATED, task_cpu(rp->task), NULL);
    kfree(rp);
}

static void sched_schedule(struct ipanema_policy *policy, unsigned int cpu)
{
    struct task_struct *next;
    struct rusty_process *rp;
  
    next = ipanema_first_task(&per_cpu(core, cpu).rq);

    if (!next)
        return;

    rp = policy_metadata(next);
    rp->state = IPANEMA_RUNNING;
    rp->rq = NULL;
    rp->slice = SLICE_NS;

    change_state(rp->task, IPANEMA_RUNNING, task_cpu(rp->task), NULL);
}


/* ------------------------ Module related functions ------------------------ */
static int fifo_init(struct ipanema_policy *policy)
{
    return 0;
}

static int fifo_free_metadata(struct ipanema_policy *policy)
{
    if (policy->data)
        kfree(policy->data);
    return 0;   
}

static int fifo_can_be_default(struct ipanema_policy *policy)
{
    return 1;
}

static bool fifo_attach(struct ipanema_policy *policy, struct task_struct *task,
                        char *command)
{
    return true;
}

/* - comment from main.rs :
        high frequency tuning decisions (100 ms)
        It identifies CPUs which are not too heavily loaded and marks them 
        so they can pull tasks from  other overloaded domains on the fly

    - how they do it :
        -  for each domain 
            - for each cpu of the domain

                - get prev_stat and cur_stat
                (stats that are in proc/stat : user / nice / system / idle / iowait /
                                                irq / softirq / stolen)

                util = **calc_util (cur, prev)**
                    - calls sub_or_zero(cur , prev) for each stat *func defined in main.rs*
                        substracts prev to cur and returns 0 if underflow happens

                    - busy_use = user + sys + nice + irq + softirq + stolen
                    - total_use = idle + busy + iowait
                    - if busy_use > 0 
                        return (busy_use / total_use ).clamp(0.0, 1.0) *clamp <=> max(val, 1.0)*
                     else 
                        return 1.0
                - gets dom_util_sum and adds util to it
                - avg_tuil += util

        - avg_util /= dom_group.weight() *weight() => the number of cpus in the domain*

        - fully_utilized = avg_util >= 0.9999 *fully_utilized is a bool of the tuner object*

        - clears direct_greedy and kick_greedy cpumasks

        for each domain 
            - calculate the domain avg_util
                if no active cpu then util = 0
                otherwise dom_util_sum / nb of active cpus
    
            - enable_direct = direct_greedy_under > 0.99999 || util < direct_greedy_under
            - enable_kick = kick_greedy_under > 0.99999 || util < kick_greedy_under
    
            *direct_greedy_under default value = 90.0 
                idle_cpus with utilization lower than this will get remote tasks directely pushed to them
    
            *kick_direct_under default value = 100.0
                idle cppus with utilization lower than this may get kicked to accelerate stealing when a task
                is queued on a saturated remote domain
    
            - if enable_direct 
                direct_greedy_mask |= dom.mask() (<=> cores of topology levels)
            
            - if enable_kick 
                kick_greedy_mask |= dom.mask()
        
*/

struct ipanema_module_routines rusty_routines = {
    /* Processes related functions */
    .new_prepare = sched_new_prepare,
    .new_end = sched_new_end,
    .terminate = sched_terminate,
    .new_place = sched_new_place,
    .tick = sched_tick,
    .yield = sched_yield,
    .block = sched_block,
    .unblock_prepare = sched_unblock_prepare,
    .unblock_place = sched_unblock_place,
    .unblock_end = sched_unblock_end,
    .schedule = sched_schedule,
    .balancing_select = sched_balancing_select,
 
    
    /* CPU cores functions */
    .get_core_state = sched_get_core_state,
    .core_entry = sched_core_entry,
    .core_exit = sched_core_exit,
    .newly_idle = sched_newly_idle,
    .enter_idle = sched_enter_idle,
    .exit_idle = sched_exit_idle,

    /* Policy related functions */
    .init = fifo_init,
    .free_metadata = fifo_free_metadata,
    .can_be_default = fifo_can_be_default,
    .attach = fifo_attach

};

/* kthread fun */
static int tuner_fn(void *arg)
{
    while(!kthread_should_stop()) {
        tuner_step();
        /* the tuning frequency can be modified by the user but the default is 100ms
        * see main.rs 
        */
        msleep(100); 
    }
    return 0;
}

/* init dom table
    rusty : 
        if no NUMA then just one node with all the cores so one domain
        otherwise build NUMA hierarchy from sysfs 
*/
static int __init init_sched(void)
{   
    printk("Hey using Rusty SaaKMed\n");

    int res, cpu, sibling, cpu_dom;
    struct rusty_core *c;

    /* cpumasks clean alloc */
    cpumask_t * domain_done = kzalloc(sizeof(cpumask_size()), GFP_KERNEL);
    all_cpumask = kzalloc(sizeof(cpumask_size()), GFP_KERNEL);
    direct_greedy_cpumask= kzalloc(sizeof(cpumask_size()), GFP_KERNEL);
    kick_greedy_cpumask= kzalloc(sizeof(cpumask_size()), GFP_KERNEL);

    cpumask_clear(domain_done);
    cpumask_clear(direct_greedy_cpumask);
    cpumask_clear(kick_greedy_cpumask);//not useful right now
    cpumask_clear(all_cpumask);

    nr_doms = 0;

    /* i made different loops but it would be better to have just one big loop that
    does everything  
    */

    /* loop that creates rusty cores for each CPUs */
    for_each_possible_cpu(cpu) {

        c = &per_cpu(core, cpu);
        c->id = cpu;
        c->state = IPANEMA_ACTIVE_CORE;
        init_ipanema_rq(&c->rq, FIFO, cpu, IPANEMA_READY, 0);

        /* rusty fields */
        c->dom_rr_cur = cpu;

        /* init the idle mask like rusty does
            we can't use the idle_cpu() kernel func so i did like that
            they'll all be idle at the begining and then everything will get into place on its own       
        */
        if (sched_smt_active()) {
            for_each_cpu(sibling, cpu_smt_mask(cpu)) {
                    cpumask_set_cpu(sibling, &idle_masks.smt);       
            }
        } 

        // if (idle_cpu(cpu)) {
        cpumask_set_cpu(cpu, &idle_masks.cpu);
        cpumask_set_cpu(cpu, all_cpumask);
    }

    /* loop that round-robins the domains
    * to do kind of like them, my dom_ctx contains the id of a domain and a *topology_level 
    * so i have the cores of that domain
    * maybe dom_ctx should only contains a cpumask but since i didn't really get how the topology_level works with NUMA
    * i left it like that    
    */
    for_each_possible_cpu(cpu) {

        struct topology_level* l = per_cpu(topology_levels, cpu);
        if(!l){
            return -ENOENT; 
        }
        
        /* no need to go over the same CPUs once their domain was done */
        if( cpumask_test_cpu(cpu, domain_done)) {
                continue;
        }

        /* the condition is likely wrong sorry*/
        if ( !(per_cpu(topology_levels, cpu)->flags & DOMAIN_NUMA)){
            
            /* init the table that round robins the domains */
            dom_data[nr_doms] = kzalloc(sizeof(struct dom_ctx), GFP_KERNEL);
            dom_data[nr_doms]->id = nr_doms;
            /* get the CPU's topology struct and put it in the table */
            dom_data[nr_doms]->topology = per_cpu(topology_levels, cpu);
            
            /* mark all those CPUs as done */
            cpumask_or(domain_done, domain_done, &dom_data[nr_doms]->topology->cores);

            /* for all the CPUs of a topo level update the dom_id */
            for_each_cpu(cpu_dom, &per_cpu(topology_levels, cpu)->cores) {
                struct rusty_core *c = &per_cpu(core, cpu_dom);
                c->dom_id = nr_doms;
            }
            nr_doms++;
        }
    }

    //kthread for tuner
    tuner_thread = kthread_run( tuner_fn, NULL, "tuner_thread");
   

    /* ipanema  stuff */
    rusty_policy = kzalloc(sizeof(struct ipanema_policy), GFP_KERNEL);
    rusty_policy->kmodule = THIS_MODULE;
    rusty_policy->routines = &rusty_routines;
    strncpy(rusty_policy->name, policy_name, MAX_POLICY_NAME_LEN);

    res = ipanema_add_policy(rusty_policy);

    if (res) {
        kfree(rusty_policy);
    }
    return res;

}
module_init(init_sched);

static void __exit exit_sched(void)
{
   int res;

   kthread_stop(tuner_thread);
    
    res = ipanema_remove_policy(rusty_policy);
    kfree(rusty_policy);
}
module_exit(exit_sched);