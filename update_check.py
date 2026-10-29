# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 José Serrano Guarinos

"""Checks GitHub for a newer Timeo release.

The release workflow publishes a GitHub release for every tag and marks
alpha/beta/rc tags as pre-releases, so /releases/latest is always the newest
stable version. Its tag is cached in a small JSON file so GitHub is asked at
most once a day.

Like cosmochlore.py, it knows nothing about Olex2's phil/param system:
Timeo.py passes in the cache path and what to do once a check finishes.
"""
import json
import os
import re
import threading
import time
import urllib.request

import HttpTools

REPO = 'Yluro/timeo'
RELEASES_URL = f'https://github.com/{REPO}/releases/latest'
_API_URL = f'https://api.github.com/repos/{REPO}/releases/latest'

CHECK_INTERVAL = 24 * 3600

# A thread outlives reload_all(), so a running check is found by name rather than a module global.
_THREAD_NAME = 'TimeoUpdateCheck'


def parse_version(text):
    """Returns (major, minor, patch) from '1.3.0' or 'v1.3.0', or None for anything else."""
    match = re.fullmatch(r'v?(\d+)\.(\d+)\.(\d+)', (text or '').strip())
    return tuple(int(g) for g in match.groups()) if match else None


def read_cache(path):
    """Returns the cached {'checked', 'ok', 'tag', 'url'} dict, or {} if it is missing or unreadable."""
    try:
        with open(path, encoding='utf-8') as f:
            cache = json.load(f)
    except (OSError, ValueError):
        return {}
    return cache if isinstance(cache, dict) else {}


def _write_cache(path, cache):
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(cache, f)
    os.replace(tmp, path)


def is_due(cache, now=None):
    """True if the cache was never written, or was written more than CHECK_INTERVAL ago."""
    now = time.time() if now is None else now
    checked = cache.get('checked')
    if not isinstance(checked, (int, float)):
        return True
    # A check time in the future means the clock was changed, so check again.
    return not 0 <= now - checked < CHECK_INTERVAL


def newer_release(current_version, cache):
    """Returns (tag, url) of the cached release if it is newer than `current_version`, else None."""
    latest = parse_version(cache.get('tag'))
    current = parse_version(current_version)
    if latest is None or current is None or latest <= current:
        return None
    return cache['tag'], cache.get('url') or RELEASES_URL


def fetch_latest(timeout=5):
    """Returns (tag, url) of the latest stable release on GitHub. Raises on any failure."""
    request = urllib.request.Request(_API_URL, headers={
        'Accept': 'application/vnd.github+json',
        'User-Agent': 'Timeo',
    })
    # HttpTools uses the proxy set in Olex2's settings, falling back to the system's.
    with HttpTools.make_url_call(request, http_timeout=timeout) as response:
        release = json.load(response)
    return release['tag_name'], release.get('html_url') or RELEASES_URL


def refresh_in_background(cache_path, on_done=None, timeout=5):
    """Fetches the latest release in a daemon thread, caches it, then calls `on_done()` from that thread."""
    if any(t.name == _THREAD_NAME and t.is_alive() for t in threading.enumerate()):
        return
    threading.Thread(target=_refresh, args=(cache_path, on_done, timeout),
                     name=_THREAD_NAME, daemon=True).start()


def _refresh(cache_path, on_done, timeout):
    cache = read_cache(cache_path)
    try:
        cache['tag'], cache['url'] = fetch_latest(timeout)
        cache['ok'] = True
    except Exception:
        # Offline, rate-limited or blocked: keep the last known tag and try again tomorrow.
        cache['ok'] = False
    cache['checked'] = time.time()
    try:
        _write_cache(cache_path, cache)
    except OSError:
        pass
    if on_done:
        on_done()
