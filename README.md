Tonight's Sky — Astronomy CLI

A Python command-line tool that tells you **what planets, the Moon, and bright stars are visible** from a given place and time.  
It prints altitudes, azimuths, compass directions, and constellations so you know where to look in the sky.

---

## Features

-  Input location by **city** (via OpenStreetMap) or **latitude/longitude**
-  Choose a specific date/time or use your **current local time**
-  Sun altitude check → tells you if it’s **daylight**, **twilight**, or **night**
-  Lists **planets + Moon** above your horizon with alt/az, compass direction, and constellation
-  Shows **bright stars (Hipparcos catalog)** with magnitude filter and constellations
-  Options to adjust **minimum altitude**, **magnitude limit**, and **maximum number of stars**
-  Works anywhere on Earth

---

## Installation

Clone the repo and install dependencies:

```bash
git clone https://github.com/yourusername/tonights-sky.git
cd tonights-sky
pip install skyfield astropy geopy pytz
