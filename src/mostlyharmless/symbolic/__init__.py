"""SymPy derivations backing the machine-checked lemmas.

Each certificate is one function returning a residual, and the tests assert that
the residual is zero.  A certificate is named, version-controlled and cited from
the manual by tag, so the derivation can present only its structural steps.

The tag list is in ROADMAP.md section 4; the evidence rungs (exact /
range-limited / numerical) are decision D10.  Prefer exact expansion and
coefficient comparison to ``simplify``, so that a passing check is a proof of the
instance rather than trust in a heuristic.
"""
