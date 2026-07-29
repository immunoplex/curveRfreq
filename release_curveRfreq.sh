#!/usr/bin/env bash
# =============================================================================
# release_curveRfreq.sh — lockstep, version-only release of curveRfreq.
#
# curveRfreq has NO functional changes this cycle. We bump it to 0.4.0 purely so
# the whole curveR stack shares a version and so the worker image can pin
# curveRfreq@v0.4.0 (see build-worker.yml CURVERFREQ_REF) for reproducible builds.
#
# Usage (from the curveRfreq repo root, in Git Bash / WSL / macOS / Linux):
#   ./release_curveRfreq.sh                 # version 0.4.0, tag v0.4.0
#   ./release_curveRfreq.sh 0.4.0 v0.4.0    # explicit version + tag
#
# Env toggles:
#   RSCRIPT=/path/to/Rscript.exe   # if Rscript isn't on PATH (common on Windows)
#   NO_DOCS=1                      # skip doc regeneration (nothing changed → fine)
#   NO_SITE=1                      # devtools::document() but skip pkgdown
#   NO_NEWS=1                      # skip the NEWS.md edit
#   NO_PUSH=1                      # do everything locally, don't push / release
#
# Because it's version-only, a re-run where DESCRIPTION is already 0.4.0 and NEWS
# already has the section will find nothing to commit — the script then offers to
# just tag the current HEAD, so you can (re)create the tag without an empty commit.
#
# Never force-pushes, never deletes tags.
# =============================================================================
set -euo pipefail

VERSION="${1:-0.4.0}"
TAG="${2:-v${VERSION}}"
PKG_EXPECTED="curveRfreq"

# ── Preflight ────────────────────────────────────────────────────────────────
command -v git >/dev/null || { echo "ERROR: git not found on PATH."; exit 1; }
git rev-parse --is-inside-work-tree >/dev/null 2>&1 \
  || { echo "ERROR: not inside a git repository. cd to the curveRfreq repo first."; exit 1; }
cd "$(git rev-parse --show-toplevel)"
[ -f DESCRIPTION ] || { echo "ERROR: no DESCRIPTION here — is this the package root?"; exit 1; }

PKG_NAME="$(awk -F': *' '/^Package:/{print $2; exit}' DESCRIPTION | tr -d '\r')"
if [ "$PKG_NAME" != "$PKG_EXPECTED" ]; then
  echo "WARNING: DESCRIPTION Package is '$PKG_NAME', expected '$PKG_EXPECTED'."
  read -r -p "Continue anyway? [y/N] " ans; [ "${ans:-N}" = "y" ] || exit 1
fi

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
echo "Repo:    $(git rev-parse --show-toplevel)"
echo "Package: $PKG_NAME"
echo "Branch:  $BRANCH"
echo "Version: $VERSION   Tag: $TAG"
echo

if git rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then
  echo "ERROR: tag $TAG already exists locally. Aborting."; exit 1
fi
if git ls-remote --exit-code --tags origin "$TAG" >/dev/null 2>&1; then
  echo "ERROR: tag $TAG already exists on origin. Aborting."; exit 1
fi

# Resolve Rscript (PATH → $RSCRIPT → common Windows locations).
resolve_rscript() {
  if [ -n "${RSCRIPT:-}" ]; then echo "$RSCRIPT"; return; fi
  if command -v Rscript >/dev/null 2>&1; then echo "Rscript"; return; fi
  local cand
  for cand in "/c/Program Files/R"/R-*/bin/x64/Rscript.exe \
              "/c/Program Files/R"/R-*/bin/Rscript.exe; do
    [ -x "$cand" ] && { echo "$cand"; return; }
  done
  echo ""
}
RS="$(resolve_rscript)"

# ── Bump DESCRIPTION Version ─────────────────────────────────────────────────
if [ -n "$RS" ] && \
   "$RS" -e 'quit(status = as.integer(!requireNamespace("desc", quietly=TRUE)))' >/dev/null 2>&1; then
  "$RS" -e "desc::desc_set_version('${VERSION}'); cat('DESCRIPTION Version ->', as.character(desc::desc_get_version()), '\n')"
else
  echo "(desc/Rscript not available — editing DESCRIPTION with awk)"
  awk -v v="$VERSION" 'BEGIN{done=0}
    /^Version:/ && !done {print "Version: " v; done=1; next} {print}' \
    DESCRIPTION > DESCRIPTION.tmp && mv DESCRIPTION.tmp DESCRIPTION
fi
DESC_VER="$(awk -F': *' '/^Version:/{print $2; exit}' DESCRIPTION | tr -d '\r')"
[ "$DESC_VER" = "$VERSION" ] || { echo "ERROR: DESCRIPTION Version is '$DESC_VER', expected '$VERSION'."; exit 1; }

# ── NEWS.md ──────────────────────────────────────────────────────────────────
NEWS_BODY="$(cat <<'EOF'
* Lockstep version bump — **no functional changes**. Released so the curveR stack
  shares a version and the worker image can pin `curveRfreq@v0.4.0` for
  reproducible builds. Rebuilt/verified against curveRcore 0.4.0.
EOF
)"

if [ "${NO_NEWS:-0}" != "1" ]; then
  DATE="$(date +%Y-%m-%d)"
  HEADER="# ${PKG_NAME} ${VERSION} (${DATE})"
  if [ -f NEWS.md ] && grep -q "^# ${PKG_NAME} ${VERSION}\b" NEWS.md; then
    echo "(NEWS.md already has a ${VERSION} section — leaving it as is)"
  else
    { printf '%s\n\n%s\n\n' "$HEADER" "$NEWS_BODY"; [ -f NEWS.md ] && cat NEWS.md; } > NEWS.md.tmp \
      && mv NEWS.md.tmp NEWS.md
    echo "Prepended ${VERSION} section to NEWS.md"
  fi
fi

# ── Regenerate documentation (optional; nothing changed, but keeps it tidy) ──
if [ "${NO_DOCS:-0}" = "1" ]; then
  echo "NO_DOCS=1 — skipping doc regeneration."
elif [ -z "$RS" ]; then
  echo "WARNING: Rscript not found; skipping doc regeneration."
  echo "  (Version-only release, so this is usually fine. Set RSCRIPT or use NO_DOCS=1"
  echo "   to silence this.)"
else
  echo "Regenerating roxygen docs (NAMESPACE + man/) ..."
  "$RS" -e 'devtools::document()'
  if [ "${NO_SITE:-0}" = "1" ]; then
    echo "NO_SITE=1 — skipping pkgdown site build."
  else
    echo "Building pkgdown site (non-fatal) ..."
    if ! "$RS" -e 'pkgdown::build_site(preview = FALSE)'; then
      echo "WARNING: pkgdown::build_site() failed (docs site only; package unaffected)."
      read -r -p "Continue the release without a rebuilt site? [y/N] " ans; [ "${ans:-N}" = "y" ] || exit 1
    fi
  fi
fi

# ── Review & confirm ─────────────────────────────────────────────────────────
git add -A
echo
echo "──────── staged for release ────────"
git --no-pager status --short
echo "────────────────────────────────────"
git --no-pager diff --cached -- DESCRIPTION NEWS.md || true
echo
read -r -p "Commit (if needed), tag ${TAG}, and push? [y/N] " ans
[ "${ans:-N}" = "y" ] || { echo "Aborted before commit. (Edits are staged but uncommitted.)"; exit 1; }

# ── Commit (or, if nothing changed, tag current HEAD) ────────────────────────
if git diff --cached --quiet; then
  echo "Nothing staged to commit — DESCRIPTION/NEWS already at ${VERSION}."
  read -r -p "Tag the current HEAD ($(git rev-parse --short HEAD)) as ${TAG}? [y/N] " ans2
  [ "${ans2:-N}" = "y" ] || { echo "Aborted; no tag created."; exit 1; }
else
  git commit -m "Release ${TAG}: lockstep version bump (no functional changes); rebuilt against curveRcore 0.4.0"
fi

git tag -a "$TAG" -m "${PKG_NAME} ${VERSION}

${NEWS_BODY}"
echo "Created annotated tag ${TAG}."

# ── Push ─────────────────────────────────────────────────────────────────────
if [ "${NO_PUSH:-0}" = "1" ]; then
  echo "NO_PUSH=1 set — skipping push. To push later:"
  echo "  git push origin ${BRANCH} && git push origin ${TAG}"
  exit 0
fi
git push origin "$BRANCH"
git push origin "$TAG"
echo "Pushed ${BRANCH} and ${TAG} to origin."

# ── GitHub Release (optional) ────────────────────────────────────────────────
if command -v gh >/dev/null && gh auth status >/dev/null 2>&1; then
  printf '%s\n' "$NEWS_BODY" > .release_notes.tmp
  gh release create "$TAG" \
    --title "${PKG_NAME} ${VERSION}" \
    --notes-file .release_notes.tmp \
    --target "$BRANCH"
  rm -f .release_notes.tmp
  echo "GitHub Release ${TAG} created."
else
  REMOTE_URL="$(git remote get-url origin 2>/dev/null || echo '')"
  echo "gh CLI not available/authenticated — the tag is pushed; create the Release here:"
  case "$REMOTE_URL" in
    git@github.com:*)     SLUG="${REMOTE_URL#git@github.com:}"; SLUG="${SLUG%.git}";;
    https://github.com/*) SLUG="${REMOTE_URL#https://github.com/}"; SLUG="${SLUG%.git}";;
    *) SLUG="";;
  esac
  [ -n "$SLUG" ] && echo "  https://github.com/${SLUG}/releases/new?tag=${TAG}"
fi
