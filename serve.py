#!/usr/bin/env python3
"""Dev static server with caching DISABLED.

`python -m http.server` applies heuristic caching, and browsers cache the
module graph imported by a `type:"module"` Web Worker especially aggressively —
so after editing e.g. setup4/solvers/*.js a normal reload keeps running the OLD
worker modules (you see stale errors like "M6 is 2D-2D only" after the code was
already fixed). This server sends `Cache-Control: no-store` on every response,
so every reload re-fetches fresh — no hard-refresh gymnastics needed.

Usage (from the repo root):
    python3 serve.py            # serves on http://localhost:8000
    python3 serve.py 8080       # custom port
"""
import sys
from http.server import HTTPServer, SimpleHTTPRequestHandler


class NoCacheHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        super().end_headers()


if __name__ == '__main__':
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    print(f'Serving (no-cache) on http://localhost:{port}  —  Ctrl+C to stop')
    try:
        HTTPServer(('', port), NoCacheHandler).serve_forever()
    except KeyboardInterrupt:
        print('\nstopped')
