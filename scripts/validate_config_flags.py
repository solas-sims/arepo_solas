#!/usr/bin/env python3
"""Validate config_flags.yaml against the Makefile and the Config.sh templates.

Checks:
  1. No derived flag (config_flags.yaml: derived_flags) is ever listed as a
     standalone, user-settable option in Template-Config.sh or any
     template_config*.yaml file — those flags are computed, never set
     directly.
  2. The Makefile's hand-written CONFIGVARS derivation blocks for each
     master switch / derived flag agree with config_flags.yaml, so the two
     can't silently drift apart.

Exits non-zero with a description of every violation found.
"""
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent


def load_spec():
    with open(ROOT / "config_flags.yaml") as f:
        spec = yaml.safe_load(f)
    return spec.get("master_switches", {}), spec.get("derived_flags", {})


def template_sh_flags(path):
    text = path.read_text()
    return set(re.findall(r"^#([A-Za-z_][A-Za-z0-9_]*)", text, re.MULTILINE))


def template_yaml_flags(path):
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    return {k.upper() for k in data}


def check_not_user_settable(derived_flags, errors):
    for sh_path in [ROOT / "Template-Config.sh"]:
        if not sh_path.exists():
            continue
        flags = template_sh_flags(sh_path)
        for name in derived_flags:
            if name in flags:
                errors.append(
                    f"{name} is a derived flag (config_flags.yaml) but appears as a "
                    f"standalone option in {sh_path.relative_to(ROOT)} — remove it, "
                    f"it must never be user-settable directly"
                )

    for yaml_path in ROOT.glob("template_config*.yaml"):
        flags = template_yaml_flags(yaml_path)
        for name in derived_flags:
            if name in flags:
                errors.append(
                    f"{name} is a derived flag (config_flags.yaml) but appears as a "
                    f"standalone key in {yaml_path.relative_to(ROOT)} — remove it, "
                    f"it must never be user-settable directly"
                )


def check_makefile_agreement(master_switches, derived_flags, errors):
    makefile = (ROOT / "Makefile").read_text()

    for name, expansion in master_switches.items():
        pattern = (
            rf"ifneq \(,\$\(filter {re.escape(name)},\$\(CONFIGVARS\)\)\)\s*\n"
            rf"CONFIGVARS \+= ([^\n]+)"
        )
        m = re.search(pattern, makefile)
        if not m:
            errors.append(
                f"master switch {name}: no matching 'CONFIGVARS +=' block found in "
                f"Makefile (spec/Makefile drift, or the block's form changed)"
            )
            continue
        makefile_expansion = m.group(1).split()
        if set(makefile_expansion) != set(expansion):
            errors.append(
                f"master switch {name}: config_flags.yaml says {expansion}, "
                f"Makefile expands to {makefile_expansion}"
            )

    for name, triggers in derived_flags.items():
        # Form 1: a Makefile variable "NAME = trigger1 trigger2 ..."
        m = re.search(rf"^{re.escape(name)}\s*=\s*([^\n]+)$", makefile, re.MULTILINE)
        if m:
            makefile_triggers = m.group(1).split()
        else:
            # Form 2: inline "ifneq (,$(filter t1 t2,$(CONFIGVARS)))\nCONFIGVARS += NAME"
            pattern = (
                rf"ifneq \(,\$\(filter ([^,]+),\$\(CONFIGVARS\)\)\)\s*\n"
                rf"CONFIGVARS \+= {re.escape(name)}\b"
            )
            m2 = re.search(pattern, makefile)
            if not m2:
                errors.append(
                    f"derived flag {name}: no matching Makefile block found "
                    f"(spec/Makefile drift, or the block's form changed)"
                )
                continue
            makefile_triggers = [t.strip("$()") for t in m2.group(1).split()]

        if set(makefile_triggers) != set(triggers):
            errors.append(
                f"derived flag {name}: config_flags.yaml says triggers={triggers}, "
                f"Makefile says {makefile_triggers}"
            )


def main():
    master_switches, derived_flags = load_spec()
    errors = []

    check_not_user_settable(derived_flags, errors)
    check_makefile_agreement(master_switches, derived_flags, errors)

    if errors:
        print("config_flags.yaml validation FAILED:\n")
        for e in errors:
            print(f"  - {e}")
        print(f"\n{len(errors)} problem(s) found.")
        return 1

    print(
        f"config_flags.yaml OK — {len(master_switches)} master switch(es), "
        f"{len(derived_flags)} derived flag(s), all consistent with the Makefile "
        f"and Config.sh templates."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
