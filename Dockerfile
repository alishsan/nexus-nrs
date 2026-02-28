# Deploy Nexus-NRS web dashboard. Build from repo root.
FROM clojure:lein-2.10.0
WORKDIR /app
COPY . .
WORKDIR /app/web-dashboard
RUN lein deps
EXPOSE 3000
# Railway sets PORT at runtime; app reads System/getenv "PORT"
CMD ["lein", "run"]
