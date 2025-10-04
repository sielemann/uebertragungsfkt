"use strict";

// ============================================================================
// GLOBAL STATE
// ============================================================================

const state = {
    // System parameters
    m: 0.2,      // Mass [kg]
    c: 1.0,      // Spring constant [N/m]
    d: 0.8,      // Damping constant [Ns/m]

    // System inputs
    A: 0.5,      // Excitation amplitude [m]
    f: 1.0,      // Excitation frequency [Hz]
    
    // Simulation state
    time: 0,
    isPlaying: false,
    offset: 0,   // Horizontal offset for visualization [px]
    
    // Simulation data (will store full trajectory)
    trajectory: [],
    maxTime: 10.0,
    dt: 0.01,   // Integration time step
};

// ============================================================================
// PHYSICS ENGINE - State Space Simulation
// ============================================================================

class SpringDamperMassSystem {
    constructor() {
        this.reset();
    }
    
    reset() {
        // Initial state: [position, velocity] relative to ground
        this.x = [0, 0];
        this.time = 0;
        this.history = {
            t: [],
            x3: [],  // Position
            x3d: [],  // Velocity
            u: [],   // Ground input
        };
    }
    
    // Ground motion (sinusoidal)
    groundPosition(t) {
        return state.A * Math.sin(2 * Math.PI * state.f * t);
    }
    
    // State space dynamics: dx/dt = f(x, u, t)
    dynamics(x, u, t) {
        const [x3, x3d] = x;
        const m = state.m;
        const c = state.c;
        const d = state.d;
        
        // State equations
        const dx3 = x3d;
        const dx3d = (c/m) * u - (c/m)*(m/d) * x3d - (c/m)*x3;
        
        return [dx3, dx3d];
    }
    
    // Runge-Kutta 4th order integration
    integrate(dt) {
        const x = this.x;
        const t = this.time;
        
        // k1 = f(t, x)
        const k1 = this.dynamics(x, this.groundPosition(t), t);
        
        // k2 = f(t + dt/2, x + dt/2 * k1)
        const x2 = [x[0] + 0.5*dt*k1[0], x[1] + 0.5*dt*k1[1]];
        const k2 = this.dynamics(x2, this.groundPosition(t + 0.5*dt), t + 0.5*dt);
        
        // k3 = f(t + dt/2, x + dt/2 * k2)
        const x3 = [x[0] + 0.5*dt*k2[0], x[1] + 0.5*dt*k2[1]];
        const k3 = this.dynamics(x3, this.groundPosition(t + 0.5*dt), t + 0.5*dt);
        
        // k4 = f(t + dt, x + dt * k3)
        const x4 = [x[0] + dt*k3[0], x[1] + dt*k3[1]];
        const k4 = this.dynamics(x4, this.groundPosition(t + dt), t + dt);
        
        // Update: x(t+dt) = x(t) + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        this.x[0] += dt/6 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        this.x[1] += dt/6 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        this.time += dt;
        
        // Store history
        this.history.t.push(this.time);
        this.history.x3.push(this.x[0]);
        this.history.x3d.push(this.x[1]);
        this.history.u.push(this.groundPosition(this.time));
    }
    
    // Simulate entire trajectory
    simulateTrajectory(tMax) {
        this.reset();
        
        // Store initial state at t=0
        this.history.t.push(0);
        this.history.x3.push(this.x[0]);
        this.history.x3d.push(this.x[1]);
        this.history.u.push(this.groundPosition(0));
        
        while (this.time < tMax) {
            this.integrate(state.dt);
        }
        return this.history;
    }
    
    // Get state at specific time (interpolation)
    getStateAt(t) {
        if (this.history.t.length === 0) return { x3: 0, x3d: 0, u: 0 };
        
        // Find closest index
        let idx = 0;
        for (let i = 0; i < this.history.t.length; i++) {
            if (this.history.t[i] >= t) {
                idx = i;
                break;
            }
        }
        
        return {
            x3: this.history.x3[idx] || 0,
            x3d: this.history.x3d[idx] || 0,
            u: this.history.u[idx] || 0,
        };
    }
    
    // Convert state to visualization positions (absolute coordinates)
    getVisualizationPositions(t) {
        const stateData = this.getStateAt(t);
        const x3 = stateData.x3;       // Mass position (absolute)
        const x3d = stateData.x3d;     // Mass velocity (absolute)
        const u = stateData.u;         // Moving attachment point position
        
        // Access parameters
        const m = state.m;
        const c = state.c;
        const d = state.d;

        // Fixed ground at reference position
        const x_ground = 0;
        
        // Moving attachment point (top of spring)
        const x1 = u;

        // Junction point (spring bottom / damper top)
        const x3dd = (c/m) * u - (c/m)*(m/d) * x3d - (c/m)*x3; // Eq. 4.21
        const x2 = x1 - x3dd*m/c; // Eq. 4.20
        
        // Mass position (already absolute)
        const x3_abs = x3;        
        
        return {
            ground: x_ground,  // Fixed ground position (reference)
            x1: x1,            // Moving attachment point (top of spring)
            x2: x2,            // Junction (spring bottom / damper top)
            x3: x3_abs,        // Mass (absolute position)
            u: u               // Moving attachment (same as x1)
        };
    }
}

const system = new SpringDamperMassSystem();

// ============================================================================
// RENDERING ENGINE - Canvas-based visualization
// ============================================================================

class SchematicRenderer {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        this.setupCanvas();
    }
    
    setupCanvas() {
        const dpr = window.devicePixelRatio || 1;
        const rect = this.canvas.getBoundingClientRect();
        
        this.canvas.width = rect.width * dpr;
        this.canvas.height = 0.5*rect.height * dpr;
        
        this.ctx.scale(dpr, dpr);
        
        this.width = rect.width;
        this.height = 0.5*rect.height;
    }
    
    render(t) {
        const ctx = this.ctx;
        const w = this.width;
        const h = this.height;
        
        // Clear
        ctx.clearRect(0, 0, w, h);
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, w, h);
        
        // Get visualization positions (all in absolute coordinates)
        const pos = system.getVisualizationPositions(t);
        
        // Scale: 1m = 200px, invert Y (screen coords go down, physics coords go up)
        const scale = 200;
        const centerX = w / 2;
        const baseY = h / 2;  // Reference point (vertically centered in canvas)
        
        // Convert physics positions to screen coordinates
        const groundScreenY = baseY - pos.ground * scale;      // Fixed ground
        const x1ScreenY = baseY - pos.x1 * scale;              // Moving attachment
        const x2ScreenY = baseY - pos.x2 * scale;              // Junction
        const x3ScreenY = baseY - pos.x3 * scale;              // Mass
        
        const massSize = 0.15 * scale;  // Mass size in pixels
        
        // Calculate horizontal offsets based on slider
        const offset = state.offset;
        const offsetX1 = 0.5 * offset;   // x1 label
        const offsetSpring = 1.0 * offset;   // Spring
        const offsetX2 = 1.5 * offset;   // x2 label
        const offsetDamper = 2.0 * offset;   // Damper
        const offsetX3Label = 2.5 * offset;   // x3 label
        const offsetMass = 3.0 * offset;   // Mass
        
        // Draw fixed ground (reference) - no offset
        this.drawGround(centerX, groundScreenY);
        
        // Draw horizontal reference lines (behind everything)
        // x1: from rigid support to spring tip
        this.drawHorizontalReferenceLine(centerX, centerX + offsetSpring, x1ScreenY);
        // x2: from spring tip to damper tip
        this.drawHorizontalReferenceLine(centerX + offsetSpring, centerX + offsetDamper, x2ScreenY);
        // x3: from damper tip to mass
        this.drawHorizontalReferenceLine(centerX + offsetDamper, centerX + offsetMass, x3ScreenY);
        
        // Draw rigid support from ground to moving attachment - no offset for dashed line
        this.drawRigidSupport(centerX, groundScreenY, x1ScreenY, offsetX1);
        
        // Draw spring from moving attachment (x1) to junction (x2)
        this.drawSpring(centerX + offsetSpring, x1ScreenY, x2ScreenY);
        
        // Draw damper from junction (x2) to mass (x3)
        this.drawDamper(centerX + offsetDamper, x2ScreenY, x3ScreenY);
        
        // Draw mass at x3
        this.drawMass(centerX + offsetMass, x3ScreenY, massSize, false);  // Don't draw label on mass
        
        // Draw x2 junction marker
        this.drawPositionMarker(centerX + offsetX2, x2ScreenY, 'x₂');
        
        // Draw x3 label (separate from mass to allow different offset)
        this.drawPositionMarker(centerX + offsetX3Label, x3ScreenY, 'x₃', '#f39c12');
    }
    
    drawHorizontalReferenceLine(x1, x2, y) {
        const ctx = this.ctx;
        
        // Draw thin black horizontal line
        ctx.strokeStyle = '#000000';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(x1, y);
        ctx.lineTo(x2, y);
        ctx.stroke();
    }
    
    drawGround(x, y) {
        const ctx = this.ctx;
        const groundWidth = 200;
        const groundHeight = 15;
        
        // Ground rectangle
        ctx.fillStyle = '#7f8c8d';
        ctx.fillRect(x - groundWidth/2 - 0.6*groundWidth, y - groundHeight/2, groundWidth, groundHeight);
        
        // Ground hatching
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        for (let i = -5; i < 5; i++) {
            const hx = x + i * 20 - 0.6*groundWidth;
            ctx.beginPath();
            ctx.moveTo(hx, y - groundHeight/2);
            ctx.lineTo(hx + 10, y + groundHeight/2);
            ctx.stroke();
        }
    }
    
    drawRigidSupport(x, y1, y2, labelOffsetX = 0) {
        const ctx = this.ctx;
        
        // Draw rigid line from ground to moving attachment - stays at centerX
        ctx.strokeStyle = '#7f8c8d';
        ctx.lineWidth = 3;
        ctx.setLineDash([8, 4]);  // Dashed line
        ctx.beginPath();
        ctx.moveTo(x, y1);
        ctx.lineTo(x, y2);
        ctx.stroke();
        ctx.setLineDash([]);  // Reset to solid
        
        // Draw larger circle at the moving attachment point (x1) - with offset
        const circleRadius = 12;
        const circleX = x + labelOffsetX;
        ctx.fillStyle = '#3498db';  // Blue to match spring
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(circleX, y2, circleRadius, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
        
        // Add "x1" label on the circle
        ctx.fillStyle = '#fff';
        ctx.font = 'bold 14px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('x₁', circleX, y2);
    }
    
    drawSpring(x, yTop, yBottom) {
        const ctx = this.ctx;
        const coils = 10;
        const amplitude = 20;
        
        ctx.strokeStyle = '#3498db';
        ctx.lineWidth = 3;
        ctx.beginPath();
        
        const length = yBottom - yTop;
        const segmentHeight = length / (coils * 2);
        
        ctx.moveTo(x, yTop);
        
        for (let i = 0; i <= coils * 2; i++) {
            const yPos = yTop + i * segmentHeight;
            const xOffset = (i % 2 === 0) ? 0 : (i % 4 === 1 ? amplitude : -amplitude);
            ctx.lineTo(x + xOffset, yPos);
        }
        
        ctx.lineTo(x, yBottom);
        ctx.stroke();
    }
    
    drawDamper(x, yTop, yBottom) {
        const ctx = this.ctx;
        const width = 40;
        const pistonWidth = 6;
        
        // Damper cylinder
        ctx.fillStyle = '#ecf0f1';
        ctx.strokeStyle = '#34495e';
        ctx.lineWidth = 2;
        ctx.fillRect(x - width/2, yTop, width, (yBottom - yTop) * 0.3);
        ctx.strokeRect(x - width/2, yTop, width, (yBottom - yTop) * 0.3);
        
        // Piston rod
        ctx.strokeStyle = '#e74c3c';
        ctx.lineWidth = pistonWidth;
        ctx.beginPath();
        ctx.moveTo(x, yTop + (yBottom - yTop) * 0.15);
        ctx.lineTo(x, yBottom);
        ctx.stroke();
    }
    
    drawMass(x, y, size, showLabel = true) {
        const ctx = this.ctx;
        
        ctx.fillStyle = '#f39c12';
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        
        ctx.fillRect(x - size/2, y - size/2, size, size);
        ctx.strokeRect(x - size/2, y - size/2, size, size);
        
        // Label (optional)
        if (showLabel) {
            ctx.fillStyle = '#fff';
            ctx.font = 'bold 16px sans-serif';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText('x₃', x, y);
        }
    }
    
    drawPositionMarker(x, y, label, color = '#e74c3c') {
        const ctx = this.ctx;
        const radius = 10;
        
        // Draw circle
        ctx.fillStyle = color;
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(x, y, radius, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
        
        // Add label
        ctx.fillStyle = '#fff';
        ctx.font = 'bold 13px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(label, x, y);
    }
    
    drawArrow(x1, y1, x2, y2, color, label) {
        const ctx = this.ctx;
        const headLength = 10;
        const angle = Math.atan2(y2 - y1, x2 - x1);
        
        ctx.strokeStyle = color;
        ctx.fillStyle = color;
        ctx.lineWidth = 2;
        
        // Line
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
        
        // Arrow head
        ctx.beginPath();
        ctx.moveTo(x2, y2);
        ctx.lineTo(x2 - headLength * Math.cos(angle - Math.PI/6), 
                   y2 - headLength * Math.sin(angle - Math.PI/6));
        ctx.lineTo(x2 - headLength * Math.cos(angle + Math.PI/6), 
                   y2 - headLength * Math.sin(angle + Math.PI/6));
        ctx.closePath();
        ctx.fill();
        
        // Label
        ctx.font = '11px sans-serif';
        ctx.fillText(label, x2 + 5, y2);
    }
}

class TrajectoryRenderer {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        this.setupCanvas();
    }
    
    setupCanvas() {
        const dpr = window.devicePixelRatio || 1;
        const rect = this.canvas.getBoundingClientRect();
        
        this.canvas.width = rect.width * dpr;
        this.canvas.height = rect.height * dpr;
        
        this.ctx.scale(dpr, dpr);
        
        this.width = rect.width;
        this.height = rect.height;
    }
    
    render(currentTime) {
        const ctx = this.ctx;
        const w = this.width;
        const h = this.height;
        
        // Clear
        ctx.clearRect(0, 0, w, h);
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, w, h);
        
        if (system.history.t.length === 0) return;
        
        // Margins and dimensions
        const margin = { top: 40, right: 40, bottom: 50, left: 60 };
        const plotWidth = w - margin.left - margin.right;
        const plotHeight = (h - margin.top - margin.bottom - 30) / 2;
        
        // Find data range
        const tMax = Math.max(...system.history.t);
        const x3Max = Math.max(...system.history.x3.map(Math.abs));
        const x3dMax = Math.max(...system.history.x3d.map(Math.abs));
        
        // Scales
        const xScale = plotWidth / tMax;
        const y3Scale = plotHeight / (x3Max * 1.2 + 0.1);
        const y3dScale = plotHeight / (x3dMax * 1.2 + 0.1);
        
        // Plot 1: Position
        this.plotTrajectory(
            ctx,
            margin.left,
            margin.top,
            plotWidth,
            plotHeight,
            system.history.t,
            system.history.x3,
            currentTime,
            xScale,
            y3Scale,
            'Position [m]',
            '#3498db'
        );
        
        // Plot 2: Velocity
        this.plotTrajectory(
            ctx,
            margin.left,
            margin.top + plotHeight + 30,
            plotWidth,
            plotHeight,
            system.history.t,
            system.history.x3d,
            currentTime,
            xScale,
            y3dScale,
            'Geschwindigkeit [m/s]',
            '#27ae60'
        );
        
        // X-axis label
        ctx.fillStyle = '#2c3e50';
        ctx.font = '14px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('Zeit [s]', w / 2, h - 10);
    }
    
    plotTrajectory(ctx, x, y, w, h, tData, yData, currentTime, xScale, yScale, label, color) {
        // Background
        ctx.fillStyle = '#f8f9fa';
        ctx.fillRect(x, y, w, h);
        
        // Grid
        ctx.strokeStyle = '#dee2e6';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 5; i++) {
            const yPos = y + i * h / 5;
            ctx.beginPath();
            ctx.moveTo(x, yPos);
            ctx.lineTo(x + w, yPos);
            ctx.stroke();
        }
        
        // Axes
        const zeroY = y + h / 2;
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(x, zeroY);
        ctx.lineTo(x + w, zeroY);
        ctx.stroke();
        
        // Full trajectory (light)
        ctx.strokeStyle = color + '40';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        for (let i = 0; i < tData.length; i++) {
            const px = x + tData[i] * xScale;
            const py = zeroY - yData[i] * yScale;
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();
        
        // Trajectory up to current time (bold)
        const currentIdx = tData.findIndex(t => t >= currentTime);
        if (currentIdx > 0) {
            ctx.strokeStyle = color;
            ctx.lineWidth = 2.5;
            ctx.beginPath();
            for (let i = 0; i <= currentIdx; i++) {
                const px = x + tData[i] * xScale;
                const py = zeroY - yData[i] * yScale;
                if (i === 0) ctx.moveTo(px, py);
                else ctx.lineTo(px, py);
            }
            ctx.stroke();
            
            // Current point
            const px = x + tData[currentIdx] * xScale;
            const py = zeroY - yData[currentIdx] * yScale;
            ctx.fillStyle = '#e74c3c';
            ctx.beginPath();
            ctx.arc(px, py, 5, 0, 2 * Math.PI);
            ctx.fill();
        }
        
        // Current time marker
        const markerX = x + currentTime * xScale;
        ctx.strokeStyle = '#e74c3c';
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(markerX, y);
        ctx.lineTo(markerX, y + h);
        ctx.stroke();
        ctx.setLineDash([]);
        
        // Label
        ctx.fillStyle = '#2c3e50';
        ctx.font = 'bold 13px sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText(label, x + 5, y + 15);
    }
}

// ============================================================================
// ANIMATION LOOP - RequestAnimationFrame
// ============================================================================

let schematicRenderer, trajectoryRenderer;
let animationFrameId = null;
let lastTimestamp = 0;

function initRenderers() {
    const schematicCanvas = document.getElementById('schematic-canvas');
    const trajectoryCanvas = document.getElementById('trajectory-canvas');
    
    schematicRenderer = new SchematicRenderer(schematicCanvas);
    trajectoryRenderer = new TrajectoryRenderer(trajectoryCanvas);
}

function animate(timestamp) {
    if (state.isPlaying) {
        // Update time
        const deltaTime = lastTimestamp ? (timestamp - lastTimestamp) / 1000 : 0;
        lastTimestamp = timestamp;
        
        state.time = Math.min(state.time + deltaTime, state.maxTime);
        
        // Update slider
        document.getElementById('time-slider').value = state.time;
        document.getElementById('time-value').textContent = state.time.toFixed(2);
        
        // Stop at end
        if (state.time >= state.maxTime) {
            state.isPlaying = false;
            updatePlayPauseButtons();
        }
    }
    
    // Render
    render();
    
    // Continue loop
    animationFrameId = requestAnimationFrame(animate);
}

function render() {
    schematicRenderer.render(state.time);
    trajectoryRenderer.render(state.time);
}

function startAnimation() {
    if (!animationFrameId) {
        lastTimestamp = 0;
        animationFrameId = requestAnimationFrame(animate);
    }
}

function stopAnimation() {
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
        animationFrameId = null;
    }
}

// ============================================================================
// EVENT HANDLERS
// ============================================================================

function updatePlayPauseButtons() {
    const playBtn = document.getElementById('play-button');
    const pauseBtn = document.getElementById('pause-button');
    
    if (state.isPlaying) {
        playBtn.style.background = '#95a5a6';
        pauseBtn.style.background = '#3498db';
    } else {
        playBtn.style.background = '#3498db';
        pauseBtn.style.background = '#95a5a6';
    }
}

function renderEquations() {
    // Check if KaTeX is loaded
    if (typeof katex === 'undefined') {
        console.warn('KaTeX not loaded yet');
        return;
    }
    
    const m = state.m;
    const c = state.c;
    const d = state.d;
    const A = state.A;
    const f = state.f;
    
    try {
        // State vector
        katex.render(
            String.raw`\mathbf{x} = \begin{bmatrix} x_3 \\ \dot{x}_3 \end{bmatrix} = \begin{bmatrix} \text{Position} \\ \text{Geschwindigkeit} \end{bmatrix}`,
            document.getElementById('state-vector-eq'),
            { displayMode: true, throwOnError: false }
        );
        
        // System matrix (with actual numeric values)
        const a21 = -(c/m);
        const a22 = -(c/d);
        const b2 = (c/m);
        
        katex.render(
            String.raw`\dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 \\ ${a21.toFixed(3)} & ${a22.toFixed(3)} \end{bmatrix} \mathbf{x} + \begin{bmatrix} 0 \\ ${b2.toFixed(3)} \end{bmatrix} u`,
            document.getElementById('system-matrix-eq'),
            { displayMode: true, throwOnError: false }
        );
        
        // Differential equation (scalar form with numeric coefficients)
        const sign22 = a22 >= 0 ? '+' : '';
        const sign21 = a21 >= 0 ? '+' : '';
        
        katex.render(
            String.raw`\ddot{x}_3 = \frac{c}{m} u - \frac{c}{d} \dot{x}_3 - \frac{c}{m} x_3 = ${b2.toFixed(2)} \cdot u ${sign22} ${a22.toFixed(2)} \dot{x}_3 ${sign21} ${a21.toFixed(2)} x_3`,
            document.getElementById('diff-eq-scalar'),
            { displayMode: true, throwOnError: false }
        );
        
        // Input function
        katex.render(
            String.raw`u(t) = A \sin(2\pi f t) = ${A.toFixed(2)} \sin(2\pi \cdot ${f.toFixed(2)} \cdot t) \text{ m}`,
            document.getElementById('input-eq'),
            { displayMode: true, throwOnError: false }
        );
        
    } catch (error) {
        console.error('Error rendering equations:', error);
    }
}

function resimulate() {
    system.simulateTrajectory(state.maxTime);
    renderEquations();  // Update equations with new parameter values
    render();
}

function setupCollapsibleSections() {
    // Add click listeners to all collapsible headers
    const collapsibleSections = document.querySelectorAll('.collapsible');
    
    collapsibleSections.forEach(section => {
        const header = section.querySelector('.collapsible-header');
        
        header.addEventListener('click', () => {
            section.classList.toggle('collapsed');
            
            // Save state to localStorage
            const sectionId = section.querySelector('h2').textContent.trim();
            const isCollapsed = section.classList.contains('collapsed');
            localStorage.setItem(`section-${sectionId}`, isCollapsed ? 'collapsed' : 'expanded');
        });
        
        // Restore state from localStorage
        const sectionId = section.querySelector('h2').textContent.trim();
        const savedState = localStorage.getItem(`section-${sectionId}`);
        if (savedState === 'collapsed') {
            section.classList.add('collapsed');
        }
    });
}

function setupEventListeners() {
    // Parameter sliders
    document.getElementById('mass-slider').addEventListener('input', (e) => {
        state.m = parseFloat(e.target.value);
        document.getElementById('mass-value').textContent = state.m.toFixed(2);
        resimulate();
    });
    
    document.getElementById('spring-slider').addEventListener('input', (e) => {
        state.c = parseFloat(e.target.value);
        document.getElementById('spring-value').textContent = state.c.toFixed(1);
        resimulate();
    });
    
    document.getElementById('damping-slider').addEventListener('input', (e) => {
        state.d = parseFloat(e.target.value);
        document.getElementById('damping-value').textContent = state.d.toFixed(1);
        resimulate();
    });
    
    document.getElementById('amplitude-slider').addEventListener('input', (e) => {
        state.A = parseFloat(e.target.value);
        document.getElementById('amplitude-value').textContent = state.A.toFixed(1);
        resimulate();
    });
    
    document.getElementById('frequency-slider').addEventListener('input', (e) => {
        state.f = parseFloat(e.target.value);
        document.getElementById('frequency-value').textContent = state.f.toFixed(1);
        resimulate();
    });
    
    // Time slider
    document.getElementById('time-slider').addEventListener('input', (e) => {
        state.time = parseFloat(e.target.value);
        document.getElementById('time-value').textContent = state.time.toFixed(2);
        render();
    });
    
    // Offset slider
    document.getElementById('offset-slider').addEventListener('input', (e) => {
        state.offset = parseFloat(e.target.value);
        document.getElementById('offset-value').textContent = state.offset.toFixed(0);
        render();
    });
    
    // Control buttons
    document.getElementById('play-button').addEventListener('click', () => {
        state.isPlaying = true;
        updatePlayPauseButtons();
    });
    
    document.getElementById('pause-button').addEventListener('click', () => {
        state.isPlaying = false;
        lastTimestamp = 0;
        updatePlayPauseButtons();
    });
    
    document.getElementById('reset-button').addEventListener('click', () => {
        state.time = 0;
        state.isPlaying = false;
        lastTimestamp = 0;
        document.getElementById('time-slider').value = 0;
        document.getElementById('time-value').textContent = '0.0';
        updatePlayPauseButtons();
        render();
    });
    
    // Window resize
    let resizeTimeout;
    window.addEventListener('resize', () => {
        clearTimeout(resizeTimeout);
        resizeTimeout = setTimeout(() => {
            schematicRenderer.setupCanvas();
            trajectoryRenderer.setupCanvas();
            render();
        }, 200);
    });
}

// ============================================================================
// INITIALIZATION
// ============================================================================

document.addEventListener('DOMContentLoaded', () => {
    initRenderers();
    setupCollapsibleSections();
    setupEventListeners();
    resimulate();
    updatePlayPauseButtons();
    startAnimation();
});
