# 🧬 Live Molecular Dynamics Visualization 🚀

> *"i just think proteins are neat"*

watch your molecular dynamics simulations in real-time - like a livestream but it's just atoms vibing 🎬✨

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![Node](https://img.shields.io/badge/node-18+-green.svg)
![Arch](https://img.shields.io/badge/btw-i%20use%20arch-1793d1.svg)

## 🌟 Overview

real-time visualization of MD simulations without crusty software from 2003. streams OpenMM simulations straight to your browser with Babylon.js rendering.

it's like if PyMOL and Minecraft had a baby 🎮

**✨ Key Features:**
- ⚡ real-time trajectory streaming (polls file for changes)
- 🎨 CPK color-coded atomic spheres
- 🔄 automatic frame detection and loading
- 🖱️ WASD + mouse look controls (actual FPS controls!)
- 📊 simple architecture, minimal dependencies
- 🚀 no external visualization servers required
- 💾 memory-optimized (4000 frame buffer limit)
- 🎮 minecraft-style controls (pointer lock + WASD)
- 🐱 tested on a thinkpad running arch btw

## 🚀 Quick Start

### 📋 Prerequisites

- 🐍 **Python 3.9+** (mamba recommended because it's faster)
- 📦 **Node.js 18+** (arch users: `pacman -S nodejs`)
- ⚗️ **OpenMM** (via conda-forge)
- 📊 **MDTraj** (for trajectory verification)

### 💻 Installation

```bash
# 1. Clone repository 📥
git clone <your-repo-url>
cd md

# 2. Create Python environment 🐍
mamba create -n live-md -c conda-forge openmm mdtraj openmmtools -y
conda activate live-md

# 3. Install Node dependencies 📦
npm install

# 4. Configure environment (optional) ⚙️
cp .env.example .env
```

### 🎬 Running

requires 2 terminal windows (tmux recommended)

**Terminal A 🧪 - OpenMM simulation:**
```bash
conda activate live-md
python python/openmm_run_pdb.py
# protein goes brrrrr
# nvidia gpu = way faster, thinkpad cpu = go make tea ☕
```

**Terminal B 💻 - Web server:**
```bash
npm run dev
# or: npm start
# serves both the viewer AND data files (no mdsrv needed!)
```

**🎮 Browser controls:**
- navigate to http://127.0.0.1:5173 🌐
- press **H** to toggle HUD (hidden by default)
- press **Y** to toggle atom inspector
- **WASD** to move camera (FPS controls) 🎮
- **QE** for up/down (minecraft style)
- **click** to enable pointer lock + mouse look 🖱️
- **JIKL;** for keyboard camera rotation 🕹️

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
         │ served via Express 🌐
         ↓
┌─────────────────┐
│ 💻 Express      │  Serves everything!
│   Server        │  → /data/* (static files) 📂
│ (port 5173)     │  → /api/config (metadata) ⚙️
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
3. 💻 **Express server** serves frontend, data files, and config via `/api/config`
4. 🔍 **Frontend** polls trajectory file size via HEAD requests
5. 📥 When file grows, frontend fetches and parses new DCD data
6. 🎨 **Babylon.js** updates atom positions and renders frames smoothly

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

## 🔧 Troubleshooting

### 🚫 Simulation won't start

```bash
# Check OpenMM installation 🔍
python -c "import openmm; print(openmm.version.version)"

# Install openmmtools if missing 📦
pip install openmmtools

# Check GPU availability 🎮
python python/check_openmm_gpu.py
```

### 📭 No frames appearing

```bash
# Check if DCD file exists 📂
ls -lh data/traj.dcd

# Verify trajectory ✅
python verify_traj.py

# Check web server is serving files 🌐
curl -I http://127.0.0.1:5173/data/traj.dcd
```

### 🖥️ Viewer shows blank screen

- 🔍 Open console (F12) and check for errors
- 📂 Verify files are accessible: http://127.0.0.1:5173/data/topology.pdb
- 🌐 Check CORS headers in Network tab
- ⚙️ Verify `/api/config` returns correct JSON

### ⏸️ Animation not playing

- 🎮 Press **H** to show HUD, check pause button state
- 🐛 Check console for DCD parser errors
- ✅ Run `python verify_traj.py` to verify frames differ
- if still broken, restart everything (classic IT solution)

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
- Current: Simple Express server serving everything

### 🧪 Testing

```bash
# Test trajectory integrity ✅
python verify_traj.py

# Test web server 💻
curl http://127.0.0.1:5173/api/config
curl -I http://127.0.0.1:5173/data/topology.pdb
curl -I http://127.0.0.1:5173/data/traj.dcd
```

## 📚 References

- 📖 [OpenMM Documentation](http://docs.openmm.org/)
- 🎨 [Babylon.js Documentation](https://doc.babylonjs.com/)
- 📄 [CHARMM DCD Specification](https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html)

## 📜 License

MIT License - see LICENSE file for details 🎉

## 🙏 Acknowledgments

- 🧪 Built with [OpenMM](https://openmm.org/)
- 🎨 Visualization powered by [Babylon.js](https://www.babylonjs.com/)
- 🧬 Test protein from PDB (1erm beta-lactamase)
- 💻 Simple Express server for file serving
- 💡 Inspired by PyMOL and ChimeraX (but in a browser with FPS controls)
- ☕ Coded at 3am (as all good projects are)
- 🐱 Emotional support provided by cats
- 🐧 arch linux btw

---

*made with 💕 and way too little caffeine*

*if you found this useful (or cursed), consider starring ⭐*
