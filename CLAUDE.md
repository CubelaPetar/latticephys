# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Rules Set

1. Ask, don't assume. If something is unclear, ask before writing a single line. Never make silent assumptions about intent, architecture, or requirements.
2. Simplest solution first. Always implement the simplest thing that could work. Do not add abstractions or flexibility that weren't explicitly requested.
3. Don't touch unrelated code. If a file or function is not directly part of the current task, do not modify it, even if you think it could be improved.
4. Flag uncertainty explicitly. If you are not confident about an approach or technical detail, say so before proceeding. Confidence without certainty causes more damage than admitting a gap.
5. Keep docs short. Simple, clear, concise, only the most important info. Split by topic or date rather than growing one file. Describe the current state, not the history of how we got there. When something is done, shorten its section instead of piling on notes.
6. Correct and clean before new features. The code must be correct, consistent and concise before extending it. Check every change against this list and report what fails:
   - **Physics:** each result is checked against an analytic formula or an independent calculation in a test.
   - **No silent errors:** inputs that would give a wrong result without notice (e.g. a wave vector that doesn't fit the supercell) raise an error.
   - **Consistent API:** new functions follow the argument order, data formats and names of existing ones. If they can't, say so before writing them.
   - **No duplication:** logic needed in two places (library, examples, tests) lives in one library function.
   - **Verify:** run `pytest`, `ruff check` and every example headless before committing.
   - Open cleanup items in `TODO.md` come before new features.

