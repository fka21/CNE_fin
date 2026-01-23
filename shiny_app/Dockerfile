FROM rocker/shiny:4.4.2

# Install required R packages
RUN install2.r --error shiny DT readr dplyr ggplot2 scales || R -e "install.packages(c('shiny','DT','readr','dplyr','ggplot2','scales'), repos='https://cloud.r-project.org')"

# Copy entire project into the image so the container is self-contained.
# Build context should be the project root; docker-compose is configured to use that.
COPY . /srv/project

# Place the Shiny app at the location shiny-server serves
RUN rm -rf /srv/shiny-server/teleost_cne || true
RUN cp -r /srv/project/shiny_app /srv/shiny-server/teleost_cne

# Ensure correct permissions
RUN chown -R shiny:shiny /srv/shiny-server/teleost_cne /srv/project

EXPOSE 3838

USER shiny
WORKDIR /srv/shiny-server/teleost_cne

# Launch shiny-server (rocker/shiny image includes shiny-server)
CMD ["/usr/bin/shiny-server", "--no-daemon"]
