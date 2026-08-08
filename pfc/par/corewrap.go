package par

import "scp/pfc/core"

// Thin non-allocating wrappers so phases.go does not need to import
// and call with slice headers that might confuse escape analysis.

func corePassPitch(sh *Shard, law core.Law) {
	core.PassPitch(sh.Lo, law, 0, sh.NOwn)
}

func corePassFieldPrecess(sh *Shard, law core.Law) {
	core.PassFieldPrecess(sh.Lo, sh.Hi, law, 0, sh.NOwn)
}

func coreApplyHopOne(sh *Shard, s int) {
	core.ApplyHopOne(sh.Hi, sh.Slots, s)
}

func corePassCommitField(sh *Shard) {
	core.PassCommitField(sh.Lo, sh.Hi, 0, sh.NOwn)
}

func coreForceOne(sh *Shard, s int, law core.Law) {
	core.ForceOne(sh.Lo, sh.Hi, sh.Slots, s, law)
}

func corePassMotion(sh *Shard, law core.Law) {
	core.PassMotion(sh.Lo, sh.Hi, law, 0, sh.NOwn)
}

func corePassClock(sh *Shard, law core.Law) {
	core.PassClock(sh.Lo, law, 0, sh.NOwn)
}
