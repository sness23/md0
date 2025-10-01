# 🧬✨ live md visualization but make it gay 🏳️‍⚧️💕

> *"i just think proteins are neat"* - me, probably high on estrogen while configuring this at 3am

watch ur molecular dynamics simulations in real-time like ur watching a twitch stream but it's just atoms vibing 🎬✨

![License](https://img.shields.io/badge/license-MIT-pink.svg)
![Python](https://img.shields.io/badge/python-3.9+-blueviolet.svg)
![Node](https://img.shields.io/badge/node-18+-hotpink.svg)
![Vibes](https://img.shields.io/badge/vibes-immaculate-ff69b4.svg)
![Arch](https://img.shields.io/badge/btw-i%20use%20arch-1793d1.svg)

## 🌟 what even is this

okay so basically i wanted to watch proteins do their little dance moves in real-time without using some crusty visualization software from 2003. so i made this cursed abomination that streams OpenMM simulations straight to your browser because why tf not 🌐

it's like if PyMOL and Minecraft had a baby and that baby was raised by Babylon.js

**✨ features (that actually work sometimes):**
- ⚡ real-time trajectory streaming (it's literally just polling a file but shhhh)
- 🎨 CPK spheres that actually look kinda cute ngl
- 🔄 auto-detects new frames (when it feels like it)
- 🖱️ WASD controls because i'm a gamer girl and mouse look is non-negotiable
- 📊 no external servers required (we're self-hosting our proteins like true linux users)
- 🚀 no bloat, no systemd, just vibes
- 💾 doesn't instantly crash ur browser (4000 frame buffer so you don't oom)
- 🎮 literally minecraft controls. i can't believe this works
- 🐱 probably works on a thinkpad (tested on my x230 running arch btw)

## 🚀 speedrun setup (any% WR attempt)

### 📋 things you need (suffer)

- 🐍 **Python 3.9+** (i use mamba because conda is slow and i have adhd)
- 📦 **Node.js 18+** (if you're on arch just `pacman -S nodejs` bestie)
- ⚗️ **OpenMM** (for the actual science part or whatever)
- 📊 **MDTraj** (to make sure we didn't fuck it up)
- 🌐 **MDsrv** (serves files, not drama)
- 🐱 **cats** (optional but recommended for debugging)

### 💻 installation (aka dependency hell)

```bash
# 1. yoink the repo 📥
git clone <your-repo-url>
cd md

# 2. python env moment 🐍
# (using mamba bc im not a masochist)
mamba create -n live-md -c conda-forge openmm mdtraj openmmtools -y
conda activate live-md

# 3. get mdsrv 🌐
pip install mdsrv

# 4. node dependencies (bracing for impact) 📦
npm install
# if this breaks just delete node_modules and try again
# it's the classic IT solution: turn it off and on again

# 5. config stuff ⚙️
cp .env.example .env
# edit this if you're fancy, otherwise defaults r fine
```

### 🎬 actually running this thing

you're gonna need like 3 terminal windows open. yes i know it's cursed. yes i use tmux. yes we exist.

**Terminal A 🧪 - the science box:**
```bash
conda activate live-md
python python/openmm_run_pdb.py
# this is where the protein goes brrrrr
# if you have a nvidia gpu this will be way faster
# if you're on a thinkpad like me, go make tea ☕
```

**Terminal B 🌐 - file server (she's serving):**
```bash
mdsrv --cfg scripts/app.cfg
# literally just serves files
# port 8080 gang
```

**Terminal C 💻 - the web thingy:**
```bash
npm run dev
# or `npm start` if you're feeling spicy
# vite goes brrrr
```

**🎮 open browser (firefox supremacy):**
- go to http://127.0.0.1:5173 🌐
- press **H** to toggle the HUD (it's hidden by default bc aesthetic)
- press **Y** for atom inspector (raycasting is black magic)
- **WASD** to move like you're playing CS:GO 🎮
- **QE** for up/down (yes it's minecraft controls, cope)
- click to enable **pointer lock** and look around with ur mouse 🖱️
- **JIKL;** if you wanna be a vim elitist about camera controls 🕹️

honestly the controls are unhinged but they work and i'm not changing them 💅

## 🏗️ Architecture

### 🔧 Components

```
┌─────────────────┐
│  🧪 OpenMM      │  Generates trajectory frames
│     Runner      │  → data/topology.pdb (once) 📝
└────────┬────────┘  → data/traj.dcd (continuous) 🎬
         │
         │ writes to filesystem 💾
         ↓
┌─────────────────┐
│   📂 data/ dir  │  Shared filesystem
│  topology.pdb   │  🧬 Structure
│  traj.dcd       │  🎞️ Trajectory
└────────┬────────┘
         │
         │ served via HTTP 🌐
         ↓
┌─────────────────┐
│  🌐 MDsrv      │  File server
│  (port 8080)    │  → /data/* (with CORS) 🔓
└────────┬────────┘
         │
         │ proxied via Express
         ↓
┌─────────────────┐
│ 💻 Express      │  Serves frontend + config
│   Server        │  → /api/config (metadata) ⚙️
└────────┬────────┘  → / (viewer) 🎨
         │
         │ HTTP 🌐
         ↓
┌─────────────────┐
│  🎮 Web         │  Babylon.js CPK viewer
│     Frontend    │  - Parses PDB/DCD in browser 🔬
└─────────────────┘  - Polls for new frames 🔄
```

### 🌊 Data Flow

1. 🧪 **OpenMM** runs molecular dynamics simulation
2. 💾 **DCDReporter** appends frames to `data/traj.dcd` every 100 steps
3. 🌐 **MDsrv** serves files from `./data` directory at port 8080
4. 💻 **Express server** serves frontend and provides MDsrv URL via `/api/config`
5. 🔍 **Frontend** polls trajectory file size via HEAD requests
6. 📥 When file grows, frontend fetches and parses new DCD data
7. 🎨 **Babylon.js** updates atom positions and renders frames smoothly

### 📄 File Formats

- 📝 **PDB** (Protein Data Bank): Text format containing topology, atom names, and initial coordinates
- 🎞️ **DCD** (CHARMM/NAMD trajectory): Binary format with 3D coordinates for each frame

## 📁 Project Structure

```
md/
├── 📂 data/                    # Generated at runtime 🏃
│   ├── topology.pdb            # Structure file (beta-lactamase 1erm) 🧬
│   └── traj.dcd                # Trajectory frames (growing) 🎬
├── 📂 public/
│   └── index.html              # Babylon.js viewer (PDB/DCD parsers, CPK renderer) 🎨
├── 📂 python/
│   ├── openmm_run_pdb.py       # MD simulation runner 🧪
│   └── check_openmm_gpu.py     # GPU availability checker ⚡
├── 📂 scripts/
│   └── app.cfg                 # MDsrv configuration 🌐
├── server.js                   # Express server (file serving + config) 💻
├── verify_traj.py              # Trajectory validation tool ✅
├── package.json                # Node dependencies 📦
├── .env.example                # Configuration template ⚙️
├── CLAUDE.md                   # AI agent instructions 🤖
├── README.md                   # This file 📖
└── README-openmm_run.md        # OpenMM runner documentation 📚
```

## ⚙️ Configuration

### 🔧 Environment Variables

Create `.env` file (or use defaults):

```bash
# Server configuration 💻
PORT=5173                           # Web server port
DATA_DIR=./data                     # Output directory
MDSRV_URL=http://127.0.0.1:8080    # MDsrv file server
```

### 🧪 OpenMM Parameters

Edit `python/openmm_run_pdb.py`:

```python
temperature = 300 * unit.kelvin     # Simulation temperature 🌡️
friction = 1.0 / unit.picosecond    # Langevin friction
timestep = 2.0 * unit.femtoseconds  # Integration timestep ⏱️
```

DCD output frequency:
```python
simulation.reporters.append(app.DCDReporter(TRAJ_PATH, 100))  # Every 100 steps 📊
```

Simulation batch size:
```python
simulation.step(5000)  # 5000 steps per batch (~10 ps) 🚀
```

## 🎮 Usage

### 🌊 Basic Workflow

1. 🧪 **Start simulation** - Begins generating frames
2. 🌐 **Start MDsrv** - Serves data files
3. 💻 **Start web server** - Serves viewer frontend
4. 🌐 **Open browser** - Loads initial topology and trajectory
5. 🔄 **Auto-polling** - Detects new frames every 60 seconds (configurable)

### 🕹️ Controls

- ⏯️ **Play/Pause** - Toggle trajectory animation (or press button)
- 🔄 **Poll interval slider** - Adjust polling frequency (1-60 seconds)
- ⚡ **Animation speed slider** - Adjust frame playback speed (10-1000ms)
- **Keyboard controls:**
  - **H** - Toggle HUD (stats, controls) 📊
  - **Y** - Toggle atom inspector (shows element/residue at crosshair) 🔍
  - **WASD** - Move camera (forward/left/back/right) 🎮
  - **QE** - Move up/down ↕️
  - **JIKL;** - Look around (camera rotation) 👀
  - **Mouse click** - Enable pointer lock for FPS-style look 🖱️

### 🧬 Using Custom Structures

See [README-openmm_run.md](./README-openmm_run.md) for detailed instructions on:
- 📥 Loading custom PDB files
- ⚗️ Choosing force fields
- 💧 Setting up explicit solvent
- 🔧 Troubleshooting structure issues

## ✅ Verification

Check trajectory integrity:

```bash
python verify_traj.py
```

Output:
```
Loading trajectory... 📥

=== Trajectory Info ===
Frames: 6165 🎬
Atoms: 22 ⚛️
Time span: 0.00 - 6164.00 ps ⏱️
Timestep: 1.00 ps

=== First frame coordinates (first 3 atoms) ===
Atom 0: (0.178, 0.135, 0.064) nm 📍
...

✅ Trajectory is valid!
```

## 🚀 Performance

### ⚡ Current Performance (CPU)

- 🧬 **System:** Beta-lactamase (1erm, ~2400 atoms, implicit solvent)
- 🖥️ **Rate:** ~100-200 steps/second on modern CPU
- 🎬 **Frame generation:** ~0.5-1.0 fps (100 steps per frame)
- 💾 **File size:** ~30 KB per 100 frames
- 🧠 **Memory:** ~1-2MB for 4000 frames (max buffer)

### 🔥 Optimization Tips

1. 🎮 **Use GPU acceleration:**
   ```python
   platform = mm.Platform.getPlatformByName('CUDA')  # or 'OpenCL'
   properties = {'CudaPrecision': 'mixed'}
   simulation = app.Simulation(topology, system, integrator, platform, properties)
   ```
   Expected speedup: 10-100x! 🚀

2. 📉 **Reduce frame output frequency:**
   ```python
   app.DCDReporter(TRAJ_PATH, 500)  # Every 500 steps
   ```

3. 📈 **Increase batch size:**
   ```python
   simulation.step(50000)  # Larger batches = better throughput
   ```

## 🌐 Browser Compatibility

- ✅ Chrome/Edge 90+ 🎯
- ✅ Firefox 88+ 🦊
- ✅ Safari 14+ 🧭
- ⚠️ Mobile browsers (limited performance) 📱

Requires:
- 🎨 WebGL 2.0
- 📦 ES6 modules
- 🌐 Fetch API
- 💾 DataView/ArrayBuffer

## ⚠️ Known Limitations

1. 🔄 **Full trajectory reload** - Currently reloads entire DCD file when detecting growth (simple but works!)
2. ✨ **Catmull-Rom interpolation** - Smooth frame-to-frame blending (4+ frames required)
3. 💾 **Limited to 4000 frames** - Browser memory constraints (auto-managed buffer)
4. 🖥️ **CPU-only by default** - Simulation uses CPU platform for compatibility (GPU available!)
5. 📂 **Single trajectory** - No support for switching between trajectories
6. 🌊 **Trajectory transition smoothing** - 1-second blend when new frames load

## 🔧 when shit breaks (a troubleshooting guide)

### 🚫 simulation won't start (skill issue)

```bash
# check if openmm is even installed lmao 🔍
python -c "import openmm; print(openmm.version.version)"

# maybe you forgot openmmtools? 📦
pip install openmmtools

# check if ur gpu works (or if you're cpu gang) 🎮
python python/check_openmm_gpu.py
# spoiler: my thinkpad doesn't have a gpu and i'm suffering
```

### 📭 no frames appearing (the void stares back)

```bash
# is the dcd file like... existing? 📂
ls -lh data/traj.dcd
# if this returns nothing, the simulation isn't running bestie

# verify you didn't fuck it up ✅
python verify_traj.py

# is mdsrv actually running or did you forget 🌐
curl -I http://127.0.0.1:8080/data/traj.dcd
# if this 404s, you forgot terminal B
```

### 🖥️ viewer is a blank void (my gender)

- 🔍 F12 and check console (become the bug)
- 📂 can you even access the files? http://127.0.0.1:8080/data/topology.pdb
- 🌐 check CORS isn't being a little bitch (network tab)
- ⚙️ `/api/config` should return json not a 404

### ⏸️ animation machine broke

- 🎮 press **H** to show HUD, click the pause button (it might already be paused)
- 🐛 check console, maybe the DCD parser is having a moment
- ✅ run `python verify_traj.py` to make sure frames are actually different
- if it's still broken, idk, restart everything. classic tech support moment

## 🛠️ Development

### ✨ Adding Features

**🎨 New visualization modes:**
- Edit `public/index.html` sphere creation in `createAtoms()`
- Modify CPK_COLORS and CPK_RADII dictionaries
- Current: CPK spheres with translucent hydrogens!

**📄 Different parsers:**
- Add parser functions alongside `parsePDB()` and `parseDCD()`
- Update `loadAll()` to use new parsers
- Current parsers use Float32Array for memory efficiency! 💾

**🌐 Server-side processing:**
- Add routes in `server.js`
- Create API endpoints for frame data
- Current: Simple Express + MDsrv proxy architecture

### 🧪 Testing

```bash
# Test trajectory integrity ✅
python verify_traj.py

# Test web server 💻
curl http://127.0.0.1:5173/api/config
curl -I http://127.0.0.1:8080/data/topology.pdb

# Test MDsrv 🌐
curl -I http://127.0.0.1:8080/data/traj.dcd
```

## 📚 References

- 📖 [OpenMM Documentation](http://docs.openmm.org/)
- 🎨 [Babylon.js Documentation](https://doc.babylonjs.com/)
- 📄 [CHARMM DCD Specification](https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html)
- 🌐 [MDsrv Documentation](https://github.com/arose/mdsrv)

## 📜 license (boring legal stuff)

MIT License - basically do whatever you want, just don't sue me 🎉

if you use this for your PhD thesis and it breaks, that's on you bestie

## 🙏 acknowledgments (who to blame)

- 🧪 built with [OpenMM](https://openmm.org/) (they're doing god's work)
- 🎨 visualization via [Babylon.js](https://www.babylonjs.com/) (webgl goes brrr)
- 🧬 test protein from PDB (1erm beta-lactamase my beloved)
- 🌐 [MDsrv](https://github.com/arose/mdsrv) for file serving (she's serving looks)
- 💡 inspired by PyMOL and ChimeraX (but like, in a browser and gayer)
- ☕ coded at 3am fueled by spite and poor life choices
- 🐱 emotional support provided by my cat (she stepped on the keyboard twice)
- 🏳️‍⚧️ trans rights are human rights
- 🐧 arch linux btw (i use arch btw, did i mention i use arch?)

---

*made with 💕 by a caffeinated transfemme who thought "what if proteins but with minecraft controls"*

*if you found this useful or cursed (or both), consider starring ⭐*

*now go simulate some proteins you beautiful disaster 🧬✨*