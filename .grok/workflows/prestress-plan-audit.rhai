let meta = #{
    name: "prestress-plan-audit",
    description: "Adversarial audit of PRESTRESS plan, seeds, and pre-registered predictions before fleet launch",
    phases: [
        #{ title: "Parallel audit", detail: "four independent read-only auditors" },
        #{ title: "Skeptic", detail: "adversarial synthesis of blockers and go/no-go" },
    ],
};

// Read-only audit of the PRESTRESS campaign plan and on-disk data before
// any Phase-2 fleet launch. Args optional; defaults to repo root.

let root = if args == () || args.root == () { "/home/d/code/scp" } else { args.root };

let audit_schema = #{
    "type": "object",
    "required": ["role", "blockers", "warnings", "ok_findings", "recommended_fixes"],
    "properties": #{
        "role": #{ "type": "string" },
        "blockers": #{
            "type": "array",
            "maxItems": 12,
            "items": #{
                "type": "object",
                "required": ["id", "claim", "evidence", "severity"],
                "properties": #{
                    "id": #{ "type": "string" },
                    "claim": #{ "type": "string" },
                    "evidence": #{ "type": "string" },
                    "severity": #{ "type": "string", "enum": ["blocker", "warning", "info"] },
                },
            },
        },
        "warnings": #{
            "type": "array",
            "maxItems": 12,
            "items": #{ "type": "string" },
        },
        "ok_findings": #{
            "type": "array",
            "maxItems": 12,
            "items": #{ "type": "string" },
        },
        "recommended_fixes": #{
            "type": "array",
            "maxItems": 10,
            "items": #{ "type": "string" },
        },
        "wave1_ready": #{ "type": "boolean" },
        "notes": #{ "type": "string" },
    },
};

let verdict_schema = #{
    "type": "object",
    "required": ["go_no_go", "summary", "must_fix_before_wave1", "must_fix_before_wave2", "plan_edits", "data_edits", "trusted_convergences", "open_cross_checks"],
    "properties": #{
        "go_no_go": #{ "type": "string", "enum": ["GO", "GO_WITH_FIXES", "NO_GO"] },
        "summary": #{ "type": "string" },
        "must_fix_before_wave1": #{ "type": "array", "items": #{ "type": "string" } },
        "must_fix_before_wave2": #{ "type": "array", "items": #{ "type": "string" } },
        "plan_edits": #{ "type": "array", "items": #{ "type": "string" } },
        "data_edits": #{ "type": "array", "items": #{ "type": "string" } },
        "trusted_convergences": #{ "type": "array", "items": #{ "type": "string" } },
        "open_cross_checks": #{ "type": "array", "items": #{ "type": "string" } },
        "recommended_wave1_order": #{ "type": "array", "items": #{ "type": "string" } },
    },
};

phase("Parallel audit");

let common_rules = "";
common_rules += "You are a read-only auditor for the v89 PRESTRESS campaign at " + root + ". ";
common_rules += "Use ONLY tools: list_dir, read_file, grep, run_terminal_command for ls/wc/rg/git show (no writes, no sims). ";
common_rules += "Do NOT consult any version before v89. ";
common_rules += "READ ACTUAL FILES — do not answer from memory. ";
common_rules += "Primary docs: v89/RESUME.md, v89/prestress/LEDGER.md, v89/MASS.md (tail ~200 lines), ";
common_rules += "v89/prestress/THEORY.md, STABILITY.md, PLASTICITY.md, FORMFIND.md, morpho/MORPHO.md, ";
common_rules += "regress/PREDICTIONS.md, candidates/summary.tsv, morpho/shortlist.json / ledger.tsv. ";
common_rules += "Seed dirs: v89/prestress/candidates/*.net, morpho/nets/*.net, foam/. ";
common_rules += "Harness: run_net.py, score_net.py; kernel: v89/cellfab.c (kappa_plast, init=net). ";
common_rules += "Be concrete: cite file paths and numbers. Severity: blocker = cannot launch that wave safely; warning = should fix or register; info = note. ";

let p_inventory = common_rules;
p_inventory += "ROLE: seed inventory and file integrity. ";
p_inventory += "Tasks: (1) Count and list all .net files in candidates/ and morpho/nets/. ";
p_inventory += "(2) Check RESUME claim of '28 kernel-ready seeds' against disk. ";
p_inventory += "(3) For Wave 1 names (c2_cube150, c2_cube150_ctrl, c8_ring12, diode-16, ring8_m3, c5_tube6, c4_ring12_chords): which exist as .net? which missing? ";
p_inventory += "(4) Is diode-16 present? Is ring16_m8 present (morpho report mentions it)? Flight-corrected P18 seeds? ";
p_inventory += "(5) Spot-check 2-3 .net files: format V/E lines, finite coords, x in [0,1], non-empty edges. ";
p_inventory += "(6) Foam foam_s20260727.tsv exists and is non-trivial size? ";
p_inventory += "(7) runs/ directory empty or absent (zero sims yet)? ";
p_inventory += "Set wave1_ready true only if every Wave1 named candidate has a seed file OR an explicit documented proxy on disk.";

let p_predictions = common_rules;
p_predictions += "ROLE: pre-registered predictions consistency. ";
p_predictions += "Tasks: (1) Extract PREDICTIONS.md bars for cube, tube, chords, load line. ";
p_predictions += "(2) Cross-check against LEDGER/RESUME/MASS for same numbers (t≈1875 band 1250-2810; skin >4700; tube ≥4600; c0=4.25e-4). ";
p_predictions += "(3) Flag contradictions between formfind summary.tsv leak scores and regression load-line (gates don't set leak). ";
p_predictions += "(4) Theory leak ratios ring6:comp12:diode and chiral-pump nu=0.64/1.05 — are scoring instruments ready in score_net.py? ";
p_predictions += "(5) PLAST bars: gates≥0.9 by t≈200, roughness halved — recorded where? ";
p_predictions += "(6) List which pre-registered claims CANNOT be scored by score_net.py as written (missing metrics). ";
p_predictions += "wave1_ready = scoring can decide the cube showdown and load-line null.";

let p_plan = common_rules;
p_plan += "ROLE: plan coherence and wave ordering. ";
p_plan += "Tasks: (1) Summarize the campaign approach (outside-sim then every avenue) and whether Wave1-4 order is sound. ";
p_plan += "(2) Check RESUME 'pending at pause' vs git HEAD (commit a1d5fc1 claims battery 21/21 + plasticity committed) — is RESUME stale? ";
p_plan += "(3) Dependencies: does Wave1 need diode emission first? spectra snap_every=250 for chiral pump? multi-seed for claims? ";
p_plan += "(4) Plasticity kappa_plast default-off; laws_V2g must NOT silently contain kappa_plast — verify. ";
p_plan += "(5) Risk that Wave1 frozen-foam runs waste wall time if formfind already says cube ceiling 0.85 — is the showdown still the right first experiment? ";
p_plan += "(6) Convergences E+D vacuum bleed, A+C parasites — does the plan actually test the right discriminator first? ";
p_plan += "(7) Propose a corrected Wave1 launch order with 2-concurrent × 4-thread constraint. ";
p_plan += "wave1_ready if plan is coherent enough to launch with at most small doc fixes.";

let p_kernel = common_rules;
p_kernel += "ROLE: kernel and harness readiness. ";
p_kernel += "Tasks: (1) Confirm cellfab binary exists and is newer/equal to cellfab.c mtime or note rebuild needed. ";
p_kernel += "(2) Confirm kappa_plast/tau_harden and init=net present; default kappa_plast=0 path. ";
p_kernel += "(3) Read run_net.py: absolute paths, snap_every default 25000 vs RESUME chiral-pump need snap_every=250 — mismatch? ";
p_kernel += "(4) score_net.py: what columns does it write; can it score t_death, late dM/dt, gate means? ";
p_kernel += "(5) laws_V2g.cfg: no kappa_plast key. ";
p_kernel += "(6) Any unfinished battery process or runs/ partial files? ";
p_kernel += "(7) Binary rebuild + smoke checklist before fleet. ";
p_kernel += "wave1_ready if harness can run a frozen-foam candidate end-to-end.";

let jobs = [
    #{
        prompt: p_inventory,
        label: "audit:inventory",
        capability_mode: "read-only",
        agent_type: "explore",
        output_schema: audit_schema,
    },
    #{
        prompt: p_predictions,
        label: "audit:predictions",
        capability_mode: "read-only",
        agent_type: "explore",
        output_schema: audit_schema,
    },
    #{
        prompt: p_plan,
        label: "audit:plan",
        capability_mode: "read-only",
        agent_type: "explore",
        output_schema: audit_schema,
    },
    #{
        prompt: p_kernel,
        label: "audit:kernel",
        capability_mode: "read-only",
        agent_type: "explore",
        output_schema: audit_schema,
    },
];

let results = parallel(jobs);

let packed = [];
let i = 0;
let labels = ["inventory", "predictions", "plan", "kernel"];
for r in results {
    let entry = #{ label: labels[i], success: false, output: () };
    if r != () && r.success {
        entry.success = true;
        entry.output = r.output;
    } else if r != () {
        entry.success = false;
        entry.output = r.output;
    }
    packed.push(entry);
    i += 1;
}

// Persist raw panel for the human and the skeptic.
let raw_json = json_encode(packed);
let raw_path = write_scratch_file("audit_panel.json", raw_json);
log("panel written " + raw_path);

phase("Skeptic");

let sk = "";
sk += "You are the adversarial skeptic for the PRESTRESS Phase-2 fleet launch decision. ";
sk += "Repo root: " + root + ". Read-only tools only. ";
sk += "You have four auditor reports (JSON below). Independently re-check any critical claim by reading the cited files — do not trust auditors blindly. ";
sk += "Especially re-verify: (a) diode-16 / ring8_m3 presence on disk, (b) score_net metrics vs pre-registrations, ";
sk += "(c) run_net snap_every vs chiral-pump spectra need, (d) whether RESUME pending-battery is stale, (e) binary/source mtime. ";
sk += "Decide go_no_go: GO (launch Wave1 as written), GO_WITH_FIXES (launch after short listed fixes), or NO_GO (do not launch). ";
sk += "Prefer GO_WITH_FIXES over NO_GO when only missing seeds/doc staleness block; NO_GO if scoring cannot decide the showdown or kernel is broken. ";
sk += "recommended_wave1_order: concrete run names and flags in launch sequence for 2 concurrent jobs. ";
sk += "Write must_fix_before_wave1 as an actionable checklist. ";
sk += "\n\nAUDITOR REPORTS JSON:\n" + raw_json;

let syn = agent(sk, #{
    label: "audit:skeptic",
    capability_mode: "read-only",
    agent_type: "explore",
    output_schema: verdict_schema,
});

if syn == () || !syn.success || syn.output == () {
    complete(#{
        status: "panel_only",
        panel_path: raw_path,
        error: "skeptic failed",
        panel: packed,
    });
}

let v = syn.output;
// Prefer the skeptic's structured object; fall back to JSON encode so
// validate_only's canned stub (not a real schema object) still completes.
let report = "# PRESTRESS plan+data audit\n\n";
report += "Panel scratch: " + raw_path + "\n\n";
report += "Skeptic payload (json):\n";
report += json_encode(v);
report += "\n\nPanel payload (json):\n";
report += raw_json;
report += "\n";

let report_path = write_scratch_file("AUDIT_REPORT.md", report);
complete(#{
    go_no_go: if type_of(v) == "map" && v.go_no_go != () { v.go_no_go } else { "UNKNOWN" },
    summary: if type_of(v) == "map" && v.summary != () { v.summary } else { "see report" },
    report_path: report_path,
    panel_path: raw_path,
    skeptic: v,
    panel: packed,
    full_report: report,
});
