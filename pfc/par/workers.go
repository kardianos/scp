package par

import "sync"

// Procedural multi-shard runtime: at FromWorld we start one goroutine
// per shard. Step never spawns; it only publishes a phase id and waits.
// There is no pool type — the loop lives here next to the barrier.

// startWorkers launches one long-lived goroutine per shard (id 0..n-1).
// Call once from FromWorld when Workers > 1.
func (c *Cluster) startWorkers() {
	c.wmu = sync.Mutex{}
	c.wcond = sync.NewCond(&c.wmu)
	c.wgen = 0
	c.wdone = 0
	c.wquit = false
	n := c.Workers
	for id := 0; id < n; id++ {
		id := id
		go c.workerLoop(id)
	}
}

// workerLoop is the entire life of shard id's goroutine.
func (c *Cluster) workerLoop(id int) {
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
		phase := c.wphase
		c.wmu.Unlock()

		c.workerPhase(id, phase)

		c.wmu.Lock()
		c.wdone++
		if c.wdone == c.Workers {
			c.wcond.Broadcast() // wake Step
		}
		c.wmu.Unlock()
	}
}

// runPhase runs phase on every shard, then returns.
// Workers==1: inline on the caller (no goroutine).
// Workers>1: wake the N shard goroutines and wait for them.
func (c *Cluster) runPhase(phase int) {
	if c.Workers <= 1 {
		c.workerPhase(0, phase)
		return
	}
	c.wmu.Lock()
	c.wphase = phase
	c.wdone = 0
	c.wgen++
	c.wcond.Broadcast()
	for c.wdone < c.Workers {
		c.wcond.Wait()
	}
	c.wmu.Unlock()
}

// stopWorkers parks all shard goroutines. Safe to call once from Stop.
func (c *Cluster) stopWorkers() {
	if c.Workers <= 1 {
		return
	}
	c.wmu.Lock()
	c.wquit = true
	c.wcond.Broadcast()
	c.wmu.Unlock()
}
