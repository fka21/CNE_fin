# Deploying the Teleost CNE Shiny App on a Local Server

This guide explains how to host the Shiny app on a local server and share it with collaborators.

## Prerequisites

- A local Linux/Unix server (Ubuntu/CentOS/Debian recommended)
- Docker and Docker Compose installed on the server
- Server IP address or domain name (e.g., `your.domain.example` or `192.168.1.100`)
- Network access configured (firewall rules to allow port 3838 or 80/443)

## Option 1: Direct Docker Compose (Simple, No HTTPS)

Best for: internal networks, quick testing

### 1. Copy the project to your server

```bash
# On your local machine:
scp -r /Users/ferenc.kagan/Documents/Projects/CNE_fin user@your-server:/home/user/teleost_cne
```

Or use `git clone` if the project is in a Git repository:
```bash
# On the server:
cd /home/user
git clone <repo-url> teleost_cne
cd teleost_cne
```

### 2. Start the app with Docker Compose

```bash
# On the server:
cd teleost_cne/shiny_app
docker compose up --build -d
```

The `-d` flag runs in background. Check status:
```bash
docker compose ps
docker compose logs -f shiny  # view logs
```

### 3. Share the link with collaborators

Share: `http://<server-ip>:3838`

Example: `http://192.168.1.100:3838`

Or if you have a local DNS name: `http://teleost-cne.local:3838`

### 4. Stop the app

```bash
cd teleost_cne/shiny_app
docker compose down
```

---

## Option 2: Docker + Nginx Reverse Proxy + HTTPS (Recommended for Production)

Best for: external access, secure, custom domain

### 1. Copy the project to server

(same as Option 1)

### 2. Install nginx and certbot

```bash
# On the server (Ubuntu/Debian):
sudo apt update
sudo apt install -y nginx certbot python3-certbot-nginx

# Or on CentOS/RHEL:
sudo yum install -y nginx certbot python3-certbot-nginx
```

### 3. Create nginx config

Create `/etc/nginx/sites-available/teleost_cne.conf`:

```nginx
server {
    listen 80;
    server_name your.domain.example;  # Replace with your domain or IP

    # Redirect HTTP to HTTPS (after certbot setup)
    # return 301 https://$server_name$request_uri;

    location / {
        proxy_pass http://127.0.0.1:3838;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_cache_bypass $http_upgrade;
    }
}
```

### 4. Enable the nginx config

```bash
sudo ln -s /etc/nginx/sites-available/teleost_cne.conf /etc/nginx/sites-enabled/
sudo nginx -t  # test config
sudo systemctl restart nginx
```

### 5. Get TLS certificate (if using a domain)

```bash
# If domain is accessible from the internet:
sudo certbot --nginx -d your.domain.example

# Follow the prompts to auto-renew
```

After certbot runs, it will automatically update nginx to redirect HTTP → HTTPS.

### 6. Start the Shiny app (if not already running)

```bash
cd /home/user/teleost_cne/shiny_app
docker compose up --build -d
```

### 7. Share the link

For external domain: `https://your.domain.example`  
For internal IP behind nginx: `http://192.168.1.100` (if DNS is set up) or `http://192.168.1.100:80`

---

## Option 3: Systemd Service (Auto-restart on reboot)

Add this to your `docker-compose.yml`:
```yaml
services:
  shiny:
    restart: unless-stopped
    ...
```

Or create a systemd service file at `/etc/systemd/system/teleost-cne.service`:

```ini
[Unit]
Description=Teleost CNE Shiny App
After=docker.service
Requires=docker.service

[Service]
Type=oneshot
RemainAfterExit=true
WorkingDirectory=/home/user/teleost_cne/shiny_app
ExecStart=/usr/bin/docker compose up -d
ExecStop=/usr/bin/docker compose down
User=user

[Install]
WantedBy=multi-user.target
```

Then:
```bash
sudo systemctl daemon-reload
sudo systemctl enable teleost-cne
sudo systemctl start teleost-cne
```

---

## Monitoring & Logs

View logs:
```bash
docker compose logs -f shiny
```

Check if port 3838 is listening:
```bash
netstat -tlnp | grep 3838
```

Test the app locally on the server:
```bash
curl http://localhost:3838
```

---

## Backups

To back up the CSV data:
```bash
# On server:
cp /home/user/teleost_cne/output/teleost_specific_cne.tsv /backup/teleost_specific_cne_$(date +%Y%m%d).tsv

# Or use a cron job (add to crontab):
0 2 * * * cp /home/user/teleost_cne/output/teleost_specific_cne.tsv /backup/teleost_specific_cne_$(date +\%Y\%m\%d).tsv
```

---

## Troubleshooting

**Port 3838 already in use:**
```bash
# Stop any existing Shiny containers:
docker ps  # find container ID
docker stop <container-id>

# Or kill the process:
lsof -i :3838
kill -9 <PID>
```

**nginx not forwarding traffic:**
- Check nginx is running: `sudo systemctl status nginx`
- Check nginx logs: `sudo tail -f /var/log/nginx/error.log`
- Test config: `sudo nginx -t`

**Data not loading in app:**
- Check working directory in Docker: `docker compose exec shiny pwd`
- Verify CSV exists: `docker compose exec shiny ls -la /srv/project/output/teleost_specific_cne.tsv`
- Check app logs: `docker compose logs -f shiny`

---

## Summary

| Option | Access | Setup Time | HTTPS | Best For |
|--------|--------|-----------|-------|----------|
| Direct Docker (Option 1) | `http://IP:3838` | 5 min | ❌ | Testing, internal only |
| Docker + nginx (Option 2) | `http://IP` or `https://domain` | 20 min | ✅ | Production, external access |
| Systemd auto-restart (Option 3) | Depends on Option 1/2 | 10 min | ❌ | Persistent uptime |

For most users, **Option 1** (direct Docker) is easiest to get started, and **Option 2** is recommended for production/external sharing.
