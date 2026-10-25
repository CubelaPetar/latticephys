# Unchecked claims in the docs

Status: ready-for-agent

## Problem

The audit of 2026-10-25 found 52 physics claims in docstrings, docs and the README; 31 had no test and 4 only an indirect one. CLAUDE.md rule 6 now says: numbers and physics statements in docs, docstrings and example output are checked by a test, or left out. Printed values (model vs. measured) are results, not claims.

## Decisions

- Analytic or symmetry facts: add a test against the formula, symmetry or limit.
- Agreement with experiment and numbers only an example produces: remove from the prose; the example keeps printing the values.
- Qualitative behaviour: test where a cheap independent check exists, otherwise reword the docstring to say what is plotted.
- When done, README can say again that every result is checked.

Line numbers in the tickets are from the audit and may have moved.
