# Deploying Nexus-NRS web dashboard on Railway

## Option A: Using the Procfile (Railway detects it)

1. Connect your **nexus-nrs** repo to Railway (root = repo root).
2. Ensure **Build Command** installs dependencies and uses the subproject:
   - **Build:** `cd web-dashboard && lein deps`
3. **Start Command** can be left empty so Railway uses the Procfile:
   - Procfile defines: `web: cd web-dashboard && lein run`
4. Set **PORT**: Railway sets `PORT` automatically; the app reads it and listens on that port.

If your stack does not include Leiningen, use Option B (Dockerfile).

## Option B: Using the Dockerfile (recommended if Leiningen is not available)

1. Connect your **nexus-nrs** repo to Railway.
2. In the service **Settings**, set **Build** to use the Dockerfile (Railway will detect the root `Dockerfile`).
3. Do **not** set a custom Root Directory; the Dockerfile is at repo root and copies the whole project so `web-dashboard` can see `../src`.
4. **Start** is `lein run` inside `web-dashboard` (from the Dockerfile `CMD`). Railway will inject `PORT`; the app uses it.

## Port

The dashboard binds to **PORT** when set (e.g. by Railway); otherwise it falls back to port 3000. No extra env vars are required for the port.

## Summary

| Method   | Build command                 | Start command              |
|----------|-------------------------------|----------------------------|
| Procfile | `cd web-dashboard && lein deps` | (from Procfile) `cd web-dashboard && lein run` |
| Dockerfile | (Docker build)             | (from CMD) `lein run` in web-dashboard |
