package parcore

import "sync"

// Procedural multi-shard runtime: FromWorld starts ONE long-lived goroutine
// per shard (id 0..n-1). Cluster.Step never spawns; it posts a phase id and
// barriers. There is no work queue and no `go` per phase — the workers sleep
// on a condition variable between phases. This is the pfc/par/workers.go
// pattern verbatim (it proved zero-alloc and good scaling): the loop lives
// here next to the barrier, not behind a pool abstraction.

func (c *Cluster[C, S]) startWorkers() {
	c.wmu = sync.Mutex{}
	c.wcond = sync.NewCond(&c.wmu)
	c.wgen = 0
	c.wdone = 0
	c.wquit = false
	for id := 0; id < c.Workers; id++ {
		id := id
		go c.workerLoop(id)
	}
}

// workerLoop is the entire life of shard id's goroutine.
func (c *Cluster[C, S]) workerLoop(id int) {
	var seen uint64
	for {
		c.wmu.Lock()
		for !c.wquit && c.wgen == seen {
			c.wcond.Wait()
		}
		if c.wquit {
			c.wmu.Unlock()
			return
		}
		seen = c.wgen
		phaseIdx := c.curPhase // the phase this generation runs
		c.wmu.Unlock()

		c.phases[phaseIdx].Run(&c.shards[id].view)

		c.wmu.Lock()
		c.wdone++
		if c.wdone == c.Workers {
			c.wcond.Broadcast() // wake Step
		}
		c.wmu.Unlock()
	}
}

// runPhase dispatches phase idx's KindParallel Run to every shard concurrently
// and returns once all have finished. Workers==1 runs inline on the caller.
// No closure is allocated — workers read c.curPhase under the lock — so Step
// stays zero-alloc.
func (c *Cluster[C, S]) runPhase(idx int) {
	if c.Workers <= 1 {
		c.phases[idx].Run(&c.shards[0].view)
		return
	}
	c.wmu.Lock()
	c.curPhase = idx
	c.wdone = 0
	c.wgen++
	c.wcond.Broadcast()
	for c.wdone < c.Workers {
		c.wcond.Wait()
	}
	c.wmu.Unlock()
}

// stopWorkers parks all shard goroutines. Safe to call once from Stop.
func (c *Cluster[C, S]) stopWorkers() {
	if c.Workers <= 1 {
		return
	}
	c.wmu.Lock()
	c.wquit = true
	c.wcond.Broadcast()
	c.wmu.Unlock()
}

// Stop parks the per-shard goroutines.
func (c *Cluster[C, S]) Stop() { c.stopWorkers() }
