# Teleost CNE browser (Shiny)

This minimal Shiny app browses `output/teleost_specific_cne.csv` and provides filters and a downloadable filtered CSV.

Requirements
- R (>= 4.0)
- R packages: `shiny`, `DT`, `readr`, `dplyr`

Run locally

From the project root folder run in R or RStudio:

```sh
# in R console
shiny::runApp('shiny_app')
```

Or open `shiny_app/app.R` in RStudio and click "Run App".

Notes
- The app expects the CSV at `output/teleost_specific_cne.csv` relative to project root.
- If the CSV is large, consider running the app on a machine with sufficient memory.
 
Docker

To run the app inside Docker (recommended for local server hosting):

1. Change to the `shiny_app` directory:

```sh
cd shiny_app
```

2. Build and start with docker-compose:

```sh
docker compose up --build
```

3. Open the app at http://localhost:3838

Notes:
- The compose file mounts the project root into the container, so the app will find `output/teleost_specific_cne.csv` as long as it exists in the project.
- The image is based on `rocker/shiny:4.4.2`; modify the `Dockerfile` if you need additional R packages.
Notes:
- The Dockerfile now copies the entire project into the image, so the container is self-contained and does not require a bind mount.
- The image is based on `rocker/shiny:4.4.2`; modify the `Dockerfile` if you need additional R packages.

Sharing the app with collaborators

- If you run the container without a reverse proxy, share the server IP and port, e.g. `http://203.0.113.45:3838` (replace with your server's public IP).
- If you configure an nginx reverse proxy and TLS for `your.domain.example` (recommended), share `https://your.domain.example`.

Security note: without authentication the app will be publicly accessible at the link you share. If the data is sensitive, add access controls (nginx basic auth or an OAuth proxy).
