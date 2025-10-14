"use strict";

// ============================================================================
// GLOBAL STATE
// ============================================================================

const state = {
    // System parameters - Two parallel branches
    m0: 0.1,     // Bottom mass [kg]
    m1: 0.2,     // Middle mass in left branch [kg]
    c0: 2.6,     // Spring constant in right branch [N/m]
    c1: 3.6,     // Top spring in left branch [N/m]
    c2: 1.8,     // Bottom spring in left branch [N/m]
    d: 0.7,      // Damping constant in right branch [Ns/m]

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
    dt: 0.02,   // Integration time step
    
    // Trajectory display flags
    isShowingInput: false,
    isShowingAnalyticTotal: true,
    isShowingAnalyticTransient: false,
    isShowingAnalyticSteady: false,
    isShowingNumerical: false,
    
    // State-space matrices (computed from parameters)
    matrices: null,  // Will store {A, B, C, D}
    
    // Complex frequency parameters for analytic solution
    s_values: null,  // Will store array of {s, weight} for exponential basis functions
    
    // Transfer function data (computed from matrices and s_values)
    transferFunctionData: null,  // Will store precomputed TF matrices and values
};

// ============================================================================
// PHYSICS ENGINE - State Space Simulation
// ============================================================================

class TwoPathSpringDamperMassSystem {
    constructor() {
        this.reset();
    }
    
    reset() {
        // Initial state: [x1, x1_dot, x0, x0_dot]
        // x1 = position of middle mass m1
        // x0 = position of bottom mass m0
        this.x = [0, 0, 0, 0];
        this.time = 0;
        this.history = {
            t: [],
            x1: [],   // Position of mass m1
            x1d: [],  // Velocity of mass m1
            x0: [],   // Position of mass m0 (output)
            x0d: [],  // Velocity of mass m0
            u: [],    // Ground input
            y_transient: [],  // Transient response
            y_steady: [],     // Steady-state response
        };
    }
    
    // Ground motion (sinusoidal)
    groundPosition(t) {
        return state.A * Math.sin(2 * Math.PI * state.f * t);
    }
    
    // Get state-space matrices A, B, C, D from physical parameters
    // State vector: x = [x1, x1_dot, x0, x0_dot]'
    // Left branch: u -> c1 -> m1 -> c2 -> m0
    // Right branch: u -> c0 -> d -> m0
    getStateSpaceMatrices() {
        const {m0, m1, c0, c1, c2, d} = state;
        
        // A matrix (4x4)
        const A = math.matrix([
            [0,              1,              0,              0],
            [-(c1+c2)/m1,    0,              c2/m1,          0],
            [0,              0,              0,              1],
            [c2/m0,          0,        -(c0+c2)/m0,     -d/m0]
        ]);
        
        // B matrix (4x1)
        const B = math.matrix([
            [0],
            [c1/m1],
            [0],
            [c0/m0]
        ]);
        
        // C matrix (1x4) - output is position of m0
        const C = math.matrix([[0, 0, 1, 0]]);
        
        // D matrix (1x1)
        const D = math.matrix([[0]]);
        
        return {A, B, C, D};
    }
    
    // Compute transfer function data for all s values
    // This precomputes matrices and transfer functions that are parameter-dependent
    // but time-independent, so they can be reused in both simulation and rendering
    computeTransferFunctionData() {
        // Get state-space matrices and s_values from state
        const {A, B, C, D} = state.matrices;
        const s_values = state.s_values;
        
        if (!s_values || !state.matrices) {
            console.warn('Cannot compute transfer function data: matrices or s_values not initialized');
            return null;
        }
        
        const stateSize = A.size()[0];
        const I = math.identity(stateSize);
        const transferFunctionData = s_values.map(({s, weight}) => {
            const sI = math.multiply(s, I);
            const sI_minus_A = math.subtract(sI, A);
            const inv_sI_A = math.inv(sI_minus_A);
            
            // Precompute (sI-A)^{-1}B (used in transient response)
            const inv_sI_A_B = math.multiply(inv_sI_A, B);
            
            // Transfer function G(s) = C(sI-A)^{-1}B + D
            const G_s = math.add(
                math.multiply(math.multiply(C, inv_sI_A), B),
                D
            );
            
            // Extract scalar from 1x1 matrix
            const G_s_scalar = math.subset(G_s, math.index(0, 0));
            
            return {
                s: s,
                weight: weight,
                sI: sI,
                sI_minus_A: sI_minus_A,
                inv_sI_A: inv_sI_A,
                inv_sI_A_B: inv_sI_A_B,
                G_s: G_s_scalar
            };
        });
        
        // Store in state for access by both simulation and rendering
        state.transferFunctionData = transferFunctionData;
        
        return transferFunctionData;
    }
    
    // State space dynamics: dx/dt = f(x, u, t)
    // State: x = [x1, x1_dot, x0, x0_dot]
    dynamics(x, u, t) {
        const [x1, x1d, x0, x0d] = x;
        const {m0, m1, c0, c1, c2, d} = state;
        
        // Equations of motion:
        // m1 * x1_dd = -c1*(x1-u) - c2*(x1-x0)
        // m0 * x0_dd = -c2*(x0-x1) - c0*(x0-u) - d*x0_dot
        
        const x1_dd = (-(c1+c2)/m1)*x1 + (c2/m1)*x0 + (c1/m1)*u;
        const x0_dd = (c2/m0)*x1 + (-(c0+c2)/m0)*x0 + (-d/m0)*x0d + (c0/m0)*u;
        
        return [x1d, x1_dd, x0d, x0_dd];
    }
    
    // Runge-Kutta 4th order integration
    integrate(dt) {
        const x = this.x;
        const t = this.time;
        
        // k1 = f(t, x)
        const k1 = this.dynamics(x, this.groundPosition(t), t);
        
        // k2 = f(t + dt/2, x + dt/2 * k1)
        const x2 = [x[0] + 0.5*dt*k1[0], x[1] + 0.5*dt*k1[1], x[2] + 0.5*dt*k1[2], x[3] + 0.5*dt*k1[3]];
        const k2 = this.dynamics(x2, this.groundPosition(t + 0.5*dt), t + 0.5*dt);
        
        // k3 = f(t + dt/2, x + dt/2 * k2)
        const x3 = [x[0] + 0.5*dt*k2[0], x[1] + 0.5*dt*k2[1], x[2] + 0.5*dt*k2[2], x[3] + 0.5*dt*k2[3]];
        const k3 = this.dynamics(x3, this.groundPosition(t + 0.5*dt), t + 0.5*dt);
        
        // k4 = f(t + dt, x + dt * k3)
        const x4 = [x[0] + dt*k3[0], x[1] + dt*k3[1], x[2] + dt*k3[2], x[3] + dt*k3[3]];
        const k4 = this.dynamics(x4, this.groundPosition(t + dt), t + dt);
        
        // Update: x(t+dt) = x(t) + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        this.x[0] += dt/6 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        this.x[1] += dt/6 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        this.x[2] += dt/6 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
        this.x[3] += dt/6 * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
        this.time += dt;
        
        // Store history
        this.history.t.push(this.time);
        this.history.x1.push(this.x[0]);
        this.history.x1d.push(this.x[1]);
        this.history.x0.push(this.x[2]);
        this.history.x0d.push(this.x[3]);
        this.history.u.push(this.groundPosition(this.time));
    }
    
    // Simulate entire trajectory
    simulateTrajectory(tMax) {
        this.reset();
        
        // Store initial state at t=0
        this.history.t.push(0);
        this.history.x1.push(this.x[0]);
        this.history.x1d.push(this.x[1]);
        this.history.x0.push(this.x[2]);
        this.history.x0d.push(this.x[3]);
        this.history.u.push(this.groundPosition(0));
        
        // Get state-space matrices from state (computed in resimulate())
        const {A, B, C, D} = state.matrices;
        const stateSize = A.size()[0];

        // Store initial state for transient response calculation
        const x0 = math.matrix(this.x.map(val => [val]));
        
        // Get precomputed transfer function data from state
        const transferFunctionData = state.transferFunctionData;
        
        if (!transferFunctionData) {
            console.error('Transfer function data not computed. Call computeTransferFunctionData() first.');
            return this.history;
        }
        
        while (this.time < tMax) {
            // Current time
            const t = this.time;

            // Numeric integration
            this.integrate(state.dt);
            
            // Compute matrix exponential e^{At} (same for all s values)
            const At = math.multiply(A, t);
            const expAt = math.expm(At);
            
            // Sum contributions from all s values (conjugate pairs)
            let y_transient_total = math.complex(0, 0);
            let y_steady_total = math.complex(0, 0);
            
            for (const {s, weight, inv_sI_A_B, G_s} of transferFunctionData) {
                // Transient response: weight * C e^{At} [x_0 - (sI-A)^{-1}B]
                // inv_sI_A_B is precomputed in transferFunctionData
                const x0_diff = math.subtract(x0, inv_sI_A_B);
                const y_trans = math.multiply(C, math.multiply(expAt, x0_diff));
                const y_trans_scalar = math.subset(y_trans, math.index(0, 0));
                const y_trans_weighted = math.multiply(weight, y_trans_scalar);
                
                // Steady-state response: weight * G(s) * e^{st}
                const est = math.exp(math.multiply(s, t));
                const y_steady_weighted = math.multiply(math.multiply(weight, G_s), est);
                
                // Accumulate
                y_transient_total = math.add(y_transient_total, y_trans_weighted);
                y_steady_total = math.add(y_steady_total, y_steady_weighted);
            }
            
            // Extract real part (amplitude A already included in weights)
            const y_transient_real = math.re(y_transient_total);
            const y_steady_real = math.re(y_steady_total);
            
            this.history.y_transient.push(y_transient_real);
            this.history.y_steady.push(y_steady_real);
        }
        return this.history;
    }
    
    // Get state at specific time (not proper interpolation, just nearest neighbor for now)
    getStateAt(t) {
        if (this.history.t.length === 0) return { x1: 0, x1d: 0, x0: 0, x0d: 0, u: 0 };
        
        // Find closest index
        let idx = 0;
        for (let i = 0; i < this.history.t.length; i++) {
            if (this.history.t[i] >= t) {
                idx = i;
                break;
            }
        }
        
        return {
            x1: this.history.x1[idx] || 0,
            x1d: this.history.x1d[idx] || 0,
            x0: this.history.x0[idx] || 0,
            x0d: this.history.x0d[idx] || 0,
            u: this.history.u[idx] || 0,
        };
    }
    
    // Convert state to visualization positions (absolute coordinates)
    getVisualizationPositions(t) {
        const stateData = this.getStateAt(t);
        const u = stateData.u;     // Top excitation
        const x1 = stateData.x1;   // Middle mass m1 position
        const x1d = stateData.x1d; // Middle mass m1 velocity
        const x0 = stateData.x0;   // Bottom mass m0 position
        const x0d = stateData.x0d; // Bottom mass m0 velocity
        
        // Get system parameters
        const {m0, m1, c0, c1, c2, d} = state;
        
        // Compute accelerations from equations of motion
        const x1dd = (-(c1+c2)/m1)*x1 + (c2/m1)*x0 + (c1/m1)*u;
        const x0dd = (c2/m0)*x1 + (-(c0+c2)/m0)*x0 + (-d/m0)*x0d + (c0/m0)*u;
        
        // LEFT PATH: Spring c1 connects directly to mass m1
        // No intermediate junction - the bottom of c1 is at x1 (mass m1 position)
        
        // RIGHT PATH: Compute junction point between c0 and d (massless junction)
        const junction2_right = u - (d / c0) * x0d;
        
        return {
            u: u,                    // Top excitation position
            x1: x1,                  // Middle mass m1 position (also bottom of c1)
            junction2_right: junction2_right, // Junction between c0 and d (right path)
            x0: x0                   // Bottom mass m0 position (output)
        };
    }
}

const system = new TwoPathSpringDamperMassSystem();

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
        const uScreenY = baseY - pos.u * scale;                       // Top excitation
        const x1ScreenY = baseY - pos.x1 * scale;                     // Middle mass m1 (also bottom of c₁)
        const junction2RightScreenY = baseY - pos.junction2_right * scale; // Right junction (c0 bottom)
        const x0ScreenY = baseY - pos.x0 * scale;                     // Bottom mass m0
        
        const massSize = 0.15 * scale;  // Mass size in pixels
        
        // Horizontal layout (two parallel paths)
        const baseLeftX = centerX - 80;   // Left branch base x-position
        const baseRightX = centerX + 80;  // Right branch base x-position
        
        // Apply offset slider value for horizontal spacing
        const offset = state.offset;
        
        // Left branch offsets (negative = move left)
        const offsetC1 = 0 * offset;      // c₁ spring: no offset
        const offsetM1 = -1 * offset;     // m₁ mass: moved left
        const offsetC2 = -2 * offset;     // c₂ spring: moved more left
        
        // Right branch offsets (positive = move right)
        const offsetC0 = 0 * offset;      // c₀ spring: no offset
        const offsetD = 1 * offset;       // d damper: moved right
        
        // Bottom mass: no offset (stays centered)
        const offsetM0 = 0;
        
        // Compute actual x-positions
        const c1X = baseLeftX + offsetC1;
        const m1X = baseLeftX + offsetM1;
        const c2X = baseLeftX + offsetC2;
        const c0X = baseRightX + offsetC0;
        const dX = baseRightX + offsetD;
        const m0X = centerX + offsetM0;
        
        // Draw top excitation point (u)
        this.drawPositionMarker(centerX, uScreenY, 'u', '#e74c3c');
        
        // ===== LEFT BRANCH: u -> c1 -> m1 -> c2 -> m0 =====
        
        // Horizontal reference line from u to c1
        this.drawHorizontalReferenceLine(centerX, c1X, uScreenY);
        
        // Draw spring c1 from u to x1 (directly to mass m1)
        this.drawSpring(c1X, uScreenY, x1ScreenY);
        ctx.fillStyle = '#2c3e50';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('c₁', c1X - 25, (uScreenY + x1ScreenY) / 2);
        
        // Horizontal reference line from c1 to m1
        this.drawHorizontalReferenceLine(c1X, m1X, x1ScreenY);
        
        // Draw mass m1 at x1
        this.drawMass(m1X, x1ScreenY, massSize * 0.8, false);
        ctx.fillStyle = '#2c3e50';
        ctx.font = 'bold 14px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('m₁', m1X, x1ScreenY);
        
        // Horizontal reference line from m1 to c2
        this.drawHorizontalReferenceLine(m1X, c2X, x1ScreenY);
        
        // Draw spring c2 from x1 to x0
        this.drawSpring(c2X, x1ScreenY, x0ScreenY);
        ctx.fillStyle = '#2c3e50';
        ctx.font = '12px sans-serif';
        ctx.fillText('c₂', c2X - 25, (x1ScreenY + x0ScreenY) / 2);
        
        // Horizontal reference line from c2 to m0
        this.drawHorizontalReferenceLine(c2X, m0X, x0ScreenY);
        
        // ===== RIGHT BRANCH: u -> c0 -> junction2 -> d -> m0 =====
        
        // Horizontal reference line from u to c0
        this.drawHorizontalReferenceLine(centerX, c0X, uScreenY);
        
        // Draw spring c0 from u to junction2
        this.drawSpring(c0X, uScreenY, junction2RightScreenY);
        ctx.fillStyle = '#2c3e50';
        ctx.font = '12px sans-serif';
        ctx.fillText('c₀', c0X + 25, (uScreenY + junction2RightScreenY) / 2);
        
        // Horizontal reference line from junction2 to d
        this.drawHorizontalReferenceLine(c0X, dX, junction2RightScreenY);
        
        // Draw damper d from junction2 to x0
        this.drawDamper(dX, junction2RightScreenY, x0ScreenY);
        ctx.fillStyle = '#2c3e50';
        ctx.font = '12px sans-serif';
        ctx.fillText('d', dX + 25, (junction2RightScreenY + x0ScreenY) / 2);
        
        // Horizontal reference line from d to m0
        this.drawHorizontalReferenceLine(dX, m0X, x0ScreenY);
        
        // ===== BOTTOM MASS m0 (where both paths meet) =====
        this.drawMass(m0X, x0ScreenY, massSize, false);
        ctx.fillStyle = '#2c3e50';
        ctx.font = 'bold 14px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('m₀', m0X, x0ScreenY);
    }
    
    drawHorizontalReferenceLine(x1, x2, y) {
        const ctx = this.ctx;
        
        // Draw thin black horizontal line
        ctx.strokeStyle = '#999';
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
        ctx.fillStyle = '#555';
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
        ctx.strokeStyle = '#555';
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
        ctx.fillStyle = '#fff';  // White for technical style
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(circleX, y2, circleRadius, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
        
        // Add "x1" label on the circle
        ctx.fillStyle = '#000';
        ctx.font = 'bold 14px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('u', circleX, y2);
    }
    
    drawSpring(x, yTop, yBottom) {
        const ctx = this.ctx;
        const coils = 10;
        const amplitude = 20;
        
        ctx.strokeStyle = '#2c3e50';
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
        ctx.fillStyle = '#e0e0e0';
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        ctx.fillRect(x - width/2, yTop, width, (yBottom - yTop) * 0.3);
        ctx.strokeRect(x - width/2, yTop, width, (yBottom - yTop) * 0.3);
        
        // Piston rod
        ctx.strokeStyle = '#555';
        ctx.lineWidth = pistonWidth;
        ctx.beginPath();
        ctx.moveTo(x, yTop + (yBottom - yTop) * 0.15);
        ctx.lineTo(x, yBottom);
        ctx.stroke();
    }
    
    drawMass(x, y, size, showLabel = true) {
        const ctx = this.ctx;
        
        ctx.fillStyle = '#e0e0e0';
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
        this.canvas.height = 0.6*rect.height * dpr;
        
        this.ctx.scale(dpr, dpr);
        
        this.width = rect.width;
        this.height = 0.6*rect.height;
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
        const plotHeight = (h - margin.top - margin.bottom - 0*30) / 1;
        
        // Find data range
        const tMax = Math.max(...system.history.t);
        
        // Compute total analytic solution
        const y_total = system.history.y_transient.map((yt, i) => 
            yt + system.history.y_steady[i]
        );
        
        // Find max values for scaling
        const x0Max = Math.max(...system.history.x0.map(Math.abs));
        const yTransMax = Math.max(...system.history.y_transient.map(Math.abs));
        const ySteadyMax = Math.max(...system.history.y_steady.map(Math.abs));
        const yTotalMax = Math.max(...y_total.map(Math.abs));
        const uMax = Math.max(...system.history.u.map(Math.abs));
        const positionMax = Math.max(x0Max, yTransMax, ySteadyMax, yTotalMax, uMax);
        
        const x0dMax = Math.max(...system.history.x0d.map(Math.abs));
        
        // Scales
        // Y-scale: add 20% padding to max value, split between positive and negative
        const xScale = plotWidth / tMax;
        const y0Scale = (plotHeight / 2) / (positionMax * 1.2 + 0.1);
        const y0dScale = (plotHeight / 2) / (x0dMax * 1.2 + 0.1);
        
        // Build curves array based on display flags
        const curves = [];
        if (state.isShowingInput) {
            curves.push({ data: system.history.u, color: '#777777', label: 'Eingang', lineWidth: 1.0 });
        }
        if (state.isShowingAnalyticTotal) {
            curves.push({ data: y_total, color: '#27ae60', label: 'Analytisch', lineWidth: 2.0 });
        }
        if (state.isShowingAnalyticTransient) {
            curves.push({ data: system.history.y_transient, color: '#e74c3c', label: 'Übergang', lineWidth: 1 });
        }
        if (state.isShowingAnalyticSteady) {
            curves.push({ data: system.history.y_steady, color: '#3498db', label: 'Stationär', lineWidth: 1 });
        }
        if (state.isShowingNumerical) {
            curves.push({ data: system.history.x0, color: '#2c3e50', label: 'Numerisch', lineWidth: 2.0, lineDash: [5, 5] });
        }
        
        // Plot 1: Position with multiple curves
        this.plotMultipleTrajectories(
            ctx,
            margin.left,
            margin.top,
            plotWidth,
            plotHeight,
            system.history.t,
            curves,
            currentTime,
            xScale,
            y0Scale,
            ''
        );
        

        
        // Axis labels
        ctx.fillStyle = '#2c3e50';
        ctx.font = '14px sans-serif';
        
        // X-axis label
        ctx.textAlign = 'center';
        ctx.fillText('Zeit [s]', w / 2, h - 10);
        
        // Y-axis label (rotated)
        ctx.save();
        ctx.translate(15, h / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.textAlign = 'center';
        ctx.fillText('Position [m]', 0, 0);
        ctx.restore();
        
        // Numeric tick labels
        ctx.font = '12px sans-serif';
        
        // X-axis ticks
        ctx.textAlign = 'left';
        ctx.fillText('0', margin.left + 2, h - margin.bottom + 20);
        ctx.textAlign = 'right';
        ctx.fillText(tMax.toFixed(1), margin.left + plotWidth - 2, h - margin.bottom + 20);
        
        // Y-axis ticks
        // The actual displayable range matches the scale with 20% padding
        const displayMax = positionMax * 1.2 + 0.1;
        ctx.textAlign = 'right';
        ctx.fillText('0', margin.left - 5, margin.top + plotHeight / 2 + 4);
        ctx.fillText(displayMax.toFixed(2), margin.left - 5, margin.top + 4);
        ctx.fillText((-displayMax).toFixed(2), margin.left - 5, margin.top + plotHeight - 4);
    }
    
    
    
    plotMultipleTrajectories(ctx, x, y, w, h, tData, curves, currentTime, xScale, yScale, label, lineDash = []) {
        // curves is an array of {data, color, label, lineWidth}
        
        // Background
        ctx.fillStyle = '#f8f9fa';
        ctx.fillRect(x, y, w, h);
        
        // Grid
        ctx.strokeStyle = '#d0d0d0';
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
        
        // Plot each curve
        for (const curve of curves) {
            const yData = curve.data;
            const color = curve.color;
            const lineWidth = curve.lineWidth || 2;
            const lineDash = curve.lineDash || [];
            
            // Full trajectory (light)
            ctx.strokeStyle = color + '40';
            ctx.lineWidth = lineWidth * 0.6;
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
                ctx.lineWidth = lineWidth;
                ctx.setLineDash(lineDash);
                ctx.beginPath();
                for (let i = 0; i <= currentIdx; i++) {
                    const px = x + tData[i] * xScale;
                    const py = zeroY - yData[i] * yScale;
                    if (i === 0) ctx.moveTo(px, py);
                    else ctx.lineTo(px, py);
                }
                ctx.stroke();
            }
        }
        
        // Current time marker
        const markerX = x + currentTime * xScale;
        ctx.strokeStyle = '#34495e';
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(markerX, y);
        ctx.lineTo(markerX, y + h);
        ctx.stroke();
        ctx.setLineDash([]);
        
        // Main label
        ctx.fillStyle = '#2c3e50';
        ctx.font = 'bold 13px sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText(label, x + 5, y + 15);
        
        // Legend for each curve (positioned on right side)
        ctx.font = '11px sans-serif';
        ctx.textAlign = 'left';
        let legendY = y + 15;
        const legendX = x + w - 80;  // 100px from right edge
        for (const curve of curves) {
            ctx.fillStyle = curve.color;
            ctx.fillRect(legendX, legendY - 8, 15, 3);
            ctx.fillStyle = '#2c3e50';
            ctx.fillText(curve.label, legendX + 20, legendY);
            legendY += 15;
        }
    }
}


class BodeRenderer {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        this.setupCanvas();
        
        // Frequency range for Bode plot (Hz)
        this.fMin = 0.01;  // 0.01 Hz
        this.fMax = 10;    // 10 Hz
        this.numPoints = 200;
        
        // Initialize empty arrays (will be computed later when matrices are ready)
        this.frequencies = [];
        this.magnitudes = [];
        this.magnitudesDB = [];
        this.phases = [];
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
    
    computeFrequencyResponse() {
        // Get state-space matrices
        const {A, B, C, D} = state.matrices;
        
        if (!A || !B || !C || !D) {
            console.warn('State-space matrices not initialized');
            return;
        }
        
        // Generate frequency points (logarithmically spaced)
        this.frequencies = [];
        this.magnitudes = [];
        this.magnitudesDB = [];
        this.phases = [];
        
        const logFMin = Math.log10(this.fMin);
        const logFMax = Math.log10(this.fMax);
        
        for (let i = 0; i < this.numPoints; i++) {
            const logF = logFMin + (logFMax - logFMin) * i / (this.numPoints - 1);
            const f = Math.pow(10, logF);
            const omega = 2 * Math.PI * f;
            
            // Compute G(jω) = C(jωI - A)^(-1)B + D
            //  - Dont reuse values from state.transferFunctionData because they assume specific s values
            //  - Here we compute s across frequency range
            const s = math.complex(0, omega);
            const stateSize = A.size()[0];
            const I = math.identity(stateSize);
            const sI = math.multiply(s, I);
            const sI_minus_A = math.subtract(sI, A);
            
            try {
                const inv_sI_A = math.inv(sI_minus_A);
                const G_s = math.add(
                    math.multiply(math.multiply(C, inv_sI_A), B),
                    D
                );
                const G_s_scalar = math.subset(G_s, math.index(0, 0));
                
                // Magnitude and phase
                const magnitude = math.abs(G_s_scalar);
                const magnitudeDB = 20 * Math.log10(magnitude);
                const phase = math.arg(G_s_scalar) * 180 / Math.PI; // Convert to degrees
                
                this.frequencies.push(f);
                this.magnitudes.push(magnitude);
                this.magnitudesDB.push(magnitudeDB);
                this.phases.push(phase);
            } catch (error) {
                console.warn(`Could not compute transfer function at f=${f}`);
            }
        }
    }
    
    render(currentFrequency) {
        const ctx = this.ctx;
        const w = this.width;
        const h = this.height;
        
        // Clear
        ctx.clearRect(0, 0, w, h);
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, w, h);
        
        if (this.frequencies.length === 0) {
            this.computeFrequencyResponse();
            if (this.frequencies.length === 0) return;
        }
        
        // Layout: Two plots stacked vertically
        const margin = { top: 50, right: 80, bottom: 60, left: 70 };
        const plotWidth = w - margin.left - margin.right;
        const plotHeight = (h - margin.top - margin.bottom - 40) / 2; // 40px gap between plots
        
        // Plot 1: Magnitude
        const magY = margin.top;
        this.renderMagnitudePlot(ctx, margin.left, magY, plotWidth, plotHeight, currentFrequency);
        
        // Plot 2: Phase
        const phaseY = magY + plotHeight + 40;
        this.renderPhasePlot(ctx, margin.left, phaseY, plotWidth, plotHeight, currentFrequency);
    }
    
    renderMagnitudePlot(ctx, x, y, w, h, currentFrequency) {
        // Background
        ctx.fillStyle = '#f8f9fa';
        ctx.fillRect(x, y, w, h);
        
        // Find magnitude range
        const magMin = Math.min(...this.magnitudes);
        const magMax = Math.max(...this.magnitudes);
        const magRange = magMax - magMin;
        
        const magDBMin = Math.min(...this.magnitudesDB);
        const magDBMax = Math.max(...this.magnitudesDB);
        const magDBRange = magDBMax - magDBMin;
        
        // Scales
        const logFMin = Math.log10(this.fMin);
        const logFMax = Math.log10(this.fMax);
        
        // Grid (logarithmic x-axis)
        ctx.strokeStyle = '#d0d0d0';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 5; i++) {
            const yPos = y + i * h / 5;
            ctx.beginPath();
            ctx.moveTo(x, yPos);
            ctx.lineTo(x + w, yPos);
            ctx.stroke();
        }
        
        // Vertical grid lines (logarithmic decades)
        for (let decade = Math.ceil(Math.log10(this.fMin)); decade <= Math.floor(Math.log10(this.fMax)); decade++) {
            const f = Math.pow(10, decade);
            const logF = Math.log10(f);
            const xPos = x + (logF - logFMin) / (logFMax - logFMin) * w;
            ctx.strokeStyle = '#c0c0c0';
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            ctx.moveTo(xPos, y);
            ctx.lineTo(xPos, y + h);
            ctx.stroke();
        }
        
        // Plot magnitude curve (linear scale on left y-axis)
        ctx.strokeStyle = '#e74c3c';
        ctx.lineWidth = 2;
        ctx.beginPath();
        for (let i = 0; i < this.frequencies.length; i++) {
            const logF = Math.log10(this.frequencies[i]);
            const px = x + (logF - logFMin) / (logFMax - logFMin) * w;
            const py = y + h - (this.magnitudes[i] - magMin) / (magRange + 0.01) * h;
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();
        
        // Current frequency marker
        if (currentFrequency >= this.fMin && currentFrequency <= this.fMax) {
            const logF = Math.log10(currentFrequency);
            const markerX = x + (logF - logFMin) / (logFMax - logFMin) * w;
            ctx.strokeStyle = '#3498db';
            ctx.lineWidth = 2;
            ctx.setLineDash([5, 5]);
            ctx.beginPath();
            ctx.moveTo(markerX, y);
            ctx.lineTo(markerX, y + h);
            ctx.stroke();
            ctx.setLineDash([]);
        }
        
        // Axes
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(x, y);
        ctx.lineTo(x, y + h);
        ctx.lineTo(x + w, y + h);
        ctx.stroke();
        
        // Labels - Title
        ctx.fillStyle = '#2c3e50';
        ctx.font = 'bold 14px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('Magnitude', x + w / 2, y - 25);
        
        // Y-axis label (left - linear)
        ctx.save();
        ctx.translate(x - 45, y + h / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('|G(jω)|', 0, 0);
        ctx.restore();
        
        // Y-axis label (right - dB)
        ctx.save();
        ctx.translate(x + w + 60, y + h / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('|G(jω)| (dB)', 0, 0);
        ctx.restore();
        
        // Y-axis ticks (left - linear)
        ctx.fillStyle = '#2c3e50';
        ctx.font = '10px sans-serif';
        ctx.textAlign = 'right';
        for (let i = 0; i <= 4; i++) {
            const val = magMin + (magMax - magMin) * i / 4;
            const yPos = y + h - i * h / 4;
            ctx.fillText(val.toFixed(2), x - 5, yPos + 3);
        }
        
        // Y-axis ticks (right - dB)
        ctx.textAlign = 'left';
        for (let i = 0; i <= 4; i++) {
            const val = magDBMin + (magDBMax - magDBMin) * i / 4;
            const yPos = y + h - i * h / 4;
            ctx.fillText(val.toFixed(1), x + w + 5, yPos + 3);
        }
        
        // X-axis label (bottom - log)
        ctx.textAlign = 'center';
        ctx.font = '12px sans-serif';
        ctx.fillText('Frequenz f (Hz)', x + w / 2, y + h + 35);
        
        // X-axis ticks (logarithmic)
        ctx.font = '10px sans-serif';
        for (let decade = Math.ceil(Math.log10(this.fMin)); decade <= Math.floor(Math.log10(this.fMax)); decade++) {
            const f = Math.pow(10, decade);
            const logF = Math.log10(f);
            const xPos = x + (logF - logFMin) / (logFMax - logFMin) * w;
            ctx.fillText(f >= 1 ? f.toFixed(0) : f.toFixed(2), xPos, y + h + 20);
        }
    }
    
    renderPhasePlot(ctx, x, y, w, h, currentFrequency) {
        // Background
        ctx.fillStyle = '#f8f9fa';
        ctx.fillRect(x, y, w, h);
        
        // Find phase range
        const phaseMin = Math.min(...this.phases);
        const phaseMax = Math.max(...this.phases);
        const phaseRange = phaseMax - phaseMin;
        
        // Scales
        const logFMin = Math.log10(this.fMin);
        const logFMax = Math.log10(this.fMax);
        
        // Grid
        ctx.strokeStyle = '#d0d0d0';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 5; i++) {
            const yPos = y + i * h / 5;
            ctx.beginPath();
            ctx.moveTo(x, yPos);
            ctx.lineTo(x + w, yPos);
            ctx.stroke();
        }
        
        // Vertical grid lines (logarithmic decades)
        for (let decade = Math.ceil(Math.log10(this.fMin)); decade <= Math.floor(Math.log10(this.fMax)); decade++) {
            const f = Math.pow(10, decade);
            const logF = Math.log10(f);
            const xPos = x + (logF - logFMin) / (logFMax - logFMin) * w;
            ctx.strokeStyle = '#c0c0c0';
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            ctx.moveTo(xPos, y);
            ctx.lineTo(xPos, y + h);
            ctx.stroke();
        }
        
        // Plot phase curve
        ctx.strokeStyle = '#9b59b6';
        ctx.lineWidth = 2;
        ctx.beginPath();
        for (let i = 0; i < this.frequencies.length; i++) {
            const logF = Math.log10(this.frequencies[i]);
            const px = x + (logF - logFMin) / (logFMax - logFMin) * w;
            const py = y + h - (this.phases[i] - phaseMin) / (phaseRange + 0.01) * h;
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();
        
        // Current frequency marker
        if (currentFrequency >= this.fMin && currentFrequency <= this.fMax) {
            const logF = Math.log10(currentFrequency);
            const markerX = x + (logF - logFMin) / (logFMax - logFMin) * w;
            ctx.strokeStyle = '#3498db';
            ctx.lineWidth = 2;
            ctx.setLineDash([5, 5]);
            ctx.beginPath();
            ctx.moveTo(markerX, y);
            ctx.lineTo(markerX, y + h);
            ctx.stroke();
            ctx.setLineDash([]);
        }
        
        // Axes
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(x, y);
        ctx.lineTo(x, y + h);
        ctx.lineTo(x + w, y + h);
        ctx.stroke();
        
        // Labels - Title
        ctx.fillStyle = '#2c3e50';
        ctx.font = 'bold 14px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('Phase', x + w / 2, y - 25);
        
        // Y-axis label
        ctx.save();
        ctx.translate(x - 45, y + h / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('∠G(jω) (°)', 0, 0);
        ctx.restore();
        
        // Y-axis ticks
        ctx.fillStyle = '#2c3e50';
        ctx.font = '10px sans-serif';
        ctx.textAlign = 'right';
        for (let i = 0; i <= 4; i++) {
            const val = phaseMin + (phaseMax - phaseMin) * i / 4;
            const yPos = y + h - i * h / 4;
            ctx.fillText(val.toFixed(1) + '°', x - 5, yPos + 3);
        }
        
        // X-axis label (bottom - log)
        ctx.textAlign = 'center';
        ctx.font = '12px sans-serif';
        ctx.fillText('Frequenz f (Hz)', x + w / 2, y + h + 35);
        
        // X-axis ticks (logarithmic)
        ctx.font = '10px sans-serif';
        for (let decade = Math.ceil(Math.log10(this.fMin)); decade <= Math.floor(Math.log10(this.fMax)); decade++) {
            const f = Math.pow(10, decade);
            const logF = Math.log10(f);
            const xPos = x + (logF - logFMin) / (logFMax - logFMin) * w;
            ctx.fillText(f >= 1 ? f.toFixed(0) : f.toFixed(2), xPos, y + h + 20);
        }
    }
}

// ============================================================================
// ANIMATION LOOP - RequestAnimationFrame
// ============================================================================

let schematicRenderer, trajectoryRenderer, bodeRenderer;
let animationFrameId = null;
let lastTimestamp = 0;

function initRenderers() {
    const schematicCanvas = document.getElementById('schematic-canvas');
    const trajectoryCanvas = document.getElementById('trajectory-canvas');
    const bodeCanvas = document.getElementById('bode-canvas');
    
    schematicRenderer = new SchematicRenderer(schematicCanvas);
    trajectoryRenderer = new TrajectoryRenderer(trajectoryCanvas);
    bodeRenderer = new BodeRenderer(bodeCanvas);
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
        
        // Update time-dependent equations
        renderEquationsTimeDependent();
        
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
    bodeRenderer.render(state.f);
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

// Helper function to format complex numbers for display
function formatComplex(c) {
    const re = c.re;
    const im = c.im;
    if (Math.abs(re) < 1e-10 && Math.abs(im) < 1e-10) return '0';
    if (Math.abs(re) < 1e-10) {
        return im >= 0 ? `${im.toFixed(2)}j` : `${im.toFixed(2)}j`;
    }
    if (Math.abs(im) < 1e-10) return `${re.toFixed(2)}`;
    const sign = im >= 0 ? '+' : '';
    return `${re.toFixed(2)}${sign}${im.toFixed(2)}j`;
}

// Render parameter-dependent equations (called when parameters change)
function renderEquationsStatic() {
    // Check if KaTeX is loaded
    if (typeof katex === 'undefined') {
        console.warn('KaTeX not loaded yet');
        return;
    }
    
    // Check if matrices are computed
    if (!state.matrices) {
        console.warn('Matrices not computed yet');
        return;
    }
    
    const A_input = state.A;
    const f = state.f;
    
    // Extract matrix values from math.js matrices (4x4 system)
    const A_mat = state.matrices.A;
    const B_mat = state.matrices.B;
    
    const a11 = math.subset(A_mat, math.index(0, 0));
    const a12 = math.subset(A_mat, math.index(0, 1));
    const a13 = math.subset(A_mat, math.index(0, 2));
    const a14 = math.subset(A_mat, math.index(0, 3));
    const a21 = math.subset(A_mat, math.index(1, 0));
    const a22 = math.subset(A_mat, math.index(1, 1));
    const a23 = math.subset(A_mat, math.index(1, 2));
    const a24 = math.subset(A_mat, math.index(1, 3));
    const a31 = math.subset(A_mat, math.index(2, 0));
    const a32 = math.subset(A_mat, math.index(2, 1));
    const a33 = math.subset(A_mat, math.index(2, 2));
    const a34 = math.subset(A_mat, math.index(2, 3));
    const a41 = math.subset(A_mat, math.index(3, 0));
    const a42 = math.subset(A_mat, math.index(3, 1));
    const a43 = math.subset(A_mat, math.index(3, 2));
    const a44 = math.subset(A_mat, math.index(3, 3));
    
    const b1 = math.subset(B_mat, math.index(0, 0));
    const b2 = math.subset(B_mat, math.index(1, 0));
    const b3 = math.subset(B_mat, math.index(2, 0));
    const b4 = math.subset(B_mat, math.index(3, 0));
    
    try {
        // State vector
        katex.render(
            String.raw`\mathbf{x} = \begin{bmatrix} x_1 \\ \dot{x}_1 \\ x_0 \\ \dot{x}_0 \end{bmatrix}`,
            document.getElementById('state-vector-eq'),
            { displayMode: true, throwOnError: false }
        );
        
        // System matrix (both symbolic and numeric forms)
        katex.render(
            String.raw`\dot{\mathbf{x}} = \begin{bmatrix} 
            0 & 1 & 0 & 0 \\ 
            -\frac{c_1+c_2}{m_1} & 0 & \frac{c_2}{m_1} & 0 \\ 
            0 & 0 & 0 & 1 \\ 
            \frac{c_2}{m_0} & 0 & -\frac{c_0+c_2}{m_0} & -\frac{d}{m_0}
            \end{bmatrix} \mathbf{x} + \begin{bmatrix} 
            0 \\ \frac{c_1}{m_1} \\ 0 \\ \frac{c_0}{m_0}
            \end{bmatrix} u = \begin{bmatrix} 
            ${a11.toFixed(2)} & ${a12.toFixed(2)} & ${a13.toFixed(2)} & ${a14.toFixed(2)} \\ 
            ${a21.toFixed(2)} & ${a22.toFixed(2)} & ${a23.toFixed(2)} & ${a24.toFixed(2)} \\ 
            ${a31.toFixed(2)} & ${a32.toFixed(2)} & ${a33.toFixed(2)} & ${a34.toFixed(2)} \\ 
            ${a41.toFixed(2)} & ${a42.toFixed(2)} & ${a43.toFixed(2)} & ${a44.toFixed(2)} 
            \end{bmatrix} \mathbf{x} + \begin{bmatrix} 
            ${b1.toFixed(2)} \\ ${b2.toFixed(2)} \\ ${b3.toFixed(2)} \\ ${b4.toFixed(2)} 
            \end{bmatrix} u`,
            document.getElementById('system-matrix-eq'),
            { displayMode: true, throwOnError: false }
        );
        
        // Input function
        katex.render(
            String.raw`u(t) = A \sin(2\pi f t) = ${A_input.toFixed(2)} \sin(2\pi \cdot ${f.toFixed(2)} \cdot t)`,
            document.getElementById('input-eq'),
            { displayMode: true, throwOnError: false }
        );
        
        // Complex parameters for transfer function (if s_values exists)
        if (state.s_values) {
            // Build list of s values
            const s_list = state.s_values.map(({s}, idx) => {
                const s_str = formatComplex(s);
                const superscript = (idx + 1) === 1 ? '_1' : '_2';
                return `s${superscript} = ${s_str}`;
            }).join(', \\; ');

            // Build list of w values 
            const w_list = state.s_values.map(({weight}, idx) => {
                const w_str = formatComplex(weight);
                const superscript = (idx + 1) === 1 ? '_1' : '_2';
                return `w${superscript} = ${w_str}`;
            }).join(', \\; ');            
            
            katex.render(
                String.raw`${s_list} \quad \quad ${w_list}`,
                document.getElementById('complex-parameters-eq'),
                { displayMode: true, throwOnError: false }
            );
            
            // Render input in terms of complex exponentials
            const omega = 2 * Math.PI * state.f;
            katex.render(
                String.raw`u(t) = w_1 \cdot e^{s_1 t} + w_2 \cdot e^{s_2 t} = A \sin(\omega t) \quad \text{mit} \quad \omega = 2\pi f = ${omega.toFixed(2)}`,
                document.getElementById('complex-input-eq'),
                { displayMode: true, throwOnError: false }
            );
        }
        
    } catch (error) {
        console.error('Error rendering static equations:', error);
    }
}

// Render time-dependent equations (called when time changes)
function renderEquationsTimeDependent() {
    // Check if KaTeX is loaded
    if (typeof katex === 'undefined') {
        return;
    }
    
    // Check if s_values and transferFunctionData are computed
    if (!state.s_values || !state.transferFunctionData) {
        return;
    }
    
    try {
        const t = state.time;
        
        // === Input equation ===
        const term1 = math.multiply(state.s_values[0].weight, math.exp(math.multiply(state.s_values[0].s, t)));
        const term2 = math.multiply(state.s_values[1].weight, math.exp(math.multiply(state.s_values[1].s, t)));
        const sum = math.add(term1, term2);
        
        const term1_str = formatComplex(term1);
        const term2_str = formatComplex(term2);
        const sum_str = sum.re.toFixed(4);  // Should be real-valued
        
        katex.render(
            String.raw`u(t) = w_1 \cdot e^{s_1 t} + w_2 \cdot e^{s_2 t} = (${term1_str}) + (${term2_str}) = ${sum_str} \quad \text{bei} \quad t = ${state.time.toFixed(2)}`,
            document.getElementById('complex-input-now-eq'),
            { displayMode: true, throwOnError: false }
        );
        
        // === Transient response equation ===
        // Get required matrices and data
        const {A, B, C, D} = state.matrices;
        
        // Get state dimension from A matrix
        const stateDim = A.size()[0];
        
        // Create zero initial state vector [0, 0, 0, ...]
        const x0 = math.zeros(stateDim, 1);
        
        // Compute matrix exponential e^{At}
        const At = math.multiply(A, t);
        const expAt = math.expm(At);
        
        // Sum contributions from both s values (conjugate pairs)
        let y_transient_total = math.complex(0, 0);
        
        // For display: get data from first s value (will show symbolic for both)
        const tfData0 = state.transferFunctionData[0];
        const tfData1 = state.transferFunctionData[1];
        
        // Compute x0_diff for both s values
        const x0_diff0 = math.subtract(x0, tfData0.inv_sI_A_B);
        const x0_diff1 = math.subtract(x0, tfData1.inv_sI_A_B);
        
        for (const {s, weight, inv_sI_A_B} of state.transferFunctionData) {
            const x0_diff = math.subtract(x0, inv_sI_A_B);
            const y_trans = math.multiply(C, math.multiply(expAt, x0_diff));
            const y_trans_scalar = math.subset(y_trans, math.index(0, 0));
            const y_trans_weighted = math.multiply(weight, y_trans_scalar);
            y_transient_total = math.add(y_transient_total, y_trans_weighted);
        }
        
        const y_transient_real = math.re(y_transient_total);
        
        // Format matrices and vectors for LaTeX (flexible for any size)
        const formatMatrixGeneral = (mat, precision = 3) => {
            const size = mat.size();
            const rows = size[0];
            const cols = size[1];
            let latex = '\\begin{bmatrix} ';
            for (let i = 0; i < rows; i++) {
                for (let j = 0; j < cols; j++) {
                    const val = math.subset(mat, math.index(i, j));
                    latex += val.toFixed(precision);
                    if (j < cols - 1) latex += ' & ';
                }
                if (i < rows - 1) latex += ' \\\\ ';
            }
            latex += ' \\end{bmatrix}';
            return latex;
        };
        
        const formatComplexVectorGeneral = (vec) => {
            const size = vec.size();
            const rows = size[0];
            let latex = '\\begin{bmatrix} ';
            for (let i = 0; i < rows; i++) {
                const val = math.subset(vec, math.index(i, 0));
                latex += formatComplex(val);
                if (i < rows - 1) latex += ' \\\\ ';
            }
            latex += ' \\end{bmatrix}';
            return latex;
        };
        
        // Format C matrix (1 x stateDim)
        const formatCMatrix = (cmat) => {
            const cols = cmat.size()[1];
            let latex = '\\begin{bmatrix} ';
            for (let j = 0; j < cols; j++) {
                const val = math.subset(cmat, math.index(0, j));
                latex += val.toFixed(0);
                if (j < cols - 1) latex += ' & ';
            }
            latex += ' \\end{bmatrix}';
            return latex;
        };
        
        const C_str = formatCMatrix(C);
        const expAt_str = formatMatrixGeneral(expAt);
        const x0_diff0_str = formatComplexVectorGeneral(x0_diff0);
        const x0_diff1_str = formatComplexVectorGeneral(x0_diff1);
        
        // Compute intermediate result: C * e^{At} * x0_diff for each s value
        const Ce_At_x0diff0 = math.multiply(C, math.multiply(expAt, x0_diff0));
        const Ce_At_x0diff0_scalar = math.subset(Ce_At_x0diff0, math.index(0, 0));
        const Ce_At_x0diff1 = math.multiply(C, math.multiply(expAt, x0_diff1));
        const Ce_At_x0diff1_scalar = math.subset(Ce_At_x0diff1, math.index(0, 0));
        
        const w0_str = formatComplex(tfData0.weight);
        const w1_str = formatComplex(tfData1.weight);
        const Ce_At_x0diff0_str = formatComplex(Ce_At_x0diff0_scalar);
        const Ce_At_x0diff1_str = formatComplex(Ce_At_x0diff1_scalar);
        
        // Weighted contributions
        const weighted0 = math.multiply(tfData0.weight, Ce_At_x0diff0_scalar);
        const weighted1 = math.multiply(tfData1.weight, Ce_At_x0diff1_scalar);
        const weighted0_str = formatComplex(weighted0);
        const weighted1_str = formatComplex(weighted1);
        
        katex.render(
            String.raw`\begin{aligned}
            y_{trans}(t) &= \sum_{i=1}^2 w_i C e^{At} \left[ x_0 - (s_i I-A)^{-1} B \right] \\
            &= w_1 \cdot \underbrace{${C_str} ${expAt_str}}_{C e^{At}} \underbrace{ ${x0_diff0_str} }_{x_0 - (s_1 I-A)^{-1} B}  + w_2 \cdot \underbrace{${C_str} ${expAt_str}}_{C e^{At}} \underbrace{ ${x0_diff1_str} }_{x_0 - (s_2 I-A)^{-1} B} \\
            &= (${w0_str}) \cdot (${Ce_At_x0diff0_str}) + (${w1_str}) \cdot (${Ce_At_x0diff1_str}) \\
            &= (${weighted0_str}) + (${weighted1_str}) = ${y_transient_real.toFixed(4)} \quad \text{bei} \quad t = ${t.toFixed(2)}
            \end{aligned}`,
            document.getElementById('transient-response-now-eq'),
            { displayMode: true, throwOnError: false }
        );
        
        // === Steady-state response equation ===
        let y_steady_total = math.complex(0, 0);
        
        for (const {s, weight, G_s} of state.transferFunctionData) {
            const est = math.exp(math.multiply(s, t));
            const y_steady_weighted = math.multiply(math.multiply(weight, G_s), est);
            y_steady_total = math.add(y_steady_total, y_steady_weighted);
        }
        
        const y_steady_real = math.re(y_steady_total);
        
        // Format matrices for steady-state equation display
        // Helper to format number (real or complex)
        const formatNumber = (val) => {
            if (typeof val === 'number') {
                return val.toFixed(2);
            }
            if (val && typeof val === 'object' && ('re' in val || 'im' in val)) {
                return formatComplex(val);
            }
            return formatComplex(math.complex(val));
        };
        
        // Helper to format complex matrix (any size)
        const formatComplexMatrixGeneral = (mat) => {
            const size = mat.size();
            const rows = size[0];
            const cols = size[1];
            let latex = '\\begin{bmatrix} ';
            for (let i = 0; i < rows; i++) {
                for (let j = 0; j < cols; j++) {
                    const val = math.subset(mat, math.index(i, j));
                    latex += formatNumber(val);
                    if (j < cols - 1) latex += ' & ';
                }
                if (i < rows - 1) latex += ' \\\\ ';
            }
            latex += ' \\end{bmatrix}';
            return latex;
        };
        
        // Helper to format complex vector (any size)
        const formatComplexVectorGeneral_steady = (vec) => {
            const size = vec.size();
            const rows = size[0];
            let latex = '\\begin{bmatrix} ';
            for (let i = 0; i < rows; i++) {
                const val = math.subset(vec, math.index(i, 0));
                latex += formatNumber(val);
                if (i < rows - 1) latex += ' \\\\ ';
            }
            latex += ' \\end{bmatrix}';
            return latex;
        };
        
        // Get matrices C, B, and D
        const C_str_steady = formatCMatrix(C);
        const B_str_steady = formatComplexVectorGeneral_steady(B);
        const D_str_steady = formatNumber(math.subset(D, math.index(0, 0)));
        
        // Get inv_sI_A matrices for both s values
        const inv_sI_A_0_str = formatComplexMatrixGeneral(tfData0.inv_sI_A);
        const inv_sI_A_1_str = formatComplexMatrixGeneral(tfData1.inv_sI_A);
        
        // Get G(s) values
        const G_s0 = tfData0.G_s;
        const G_s1 = tfData1.G_s;
        const G_s0_str = formatComplex(G_s0);
        const G_s1_str = formatComplex(G_s1);
        
        // Compute e^{st}
        const est0 = math.exp(math.multiply(tfData0.s, t));
        const est1 = math.exp(math.multiply(tfData1.s, t));
        const est0_str = formatComplex(est0);
        const est1_str = formatComplex(est1);
        
        // Compute G(s) * e^{st}
        const G_est0 = math.multiply(G_s0, est0);
        const G_est1 = math.multiply(G_s1, est1);
        const G_est0_str = formatComplex(G_est0);
        const G_est1_str = formatComplex(G_est1);
        
        // Weighted contributions
        const steady_weighted0 = math.multiply(tfData0.weight, G_est0);
        const steady_weighted1 = math.multiply(tfData1.weight, G_est1);
        const steady_weighted0_str = formatComplex(steady_weighted0);
        const steady_weighted1_str = formatComplex(steady_weighted1);
        
        katex.render(
            String.raw`\begin{aligned}
            y_{stat}(t) &= \sum_{i=1}^2 w_i \left[ C (s_i I-A)^{-1} B + D \right] e^{s_i t} = \sum_{i=1}^2 w_i G(s_i) e^{s_i t} \\
            &= w_1 \cdot \left[ \underbrace{${C_str_steady} ${inv_sI_A_0_str} ${B_str_steady}}_{C (s_1 I-A)^{-1} B} + ${D_str_steady} \right] \cdot \underbrace{(${est0_str})}_{e^{s_1 t}} + \\
            &  \quad \quad \quad + w_2 \cdot \left[ \underbrace{${C_str_steady} ${inv_sI_A_1_str} ${B_str_steady}}_{C (s_2 I-A)^{-1} B} + ${D_str_steady} \right] \cdot \underbrace{(${est1_str})}_{e^{s_2 t}} \\
            &= w_1 \cdot (${G_s0_str}) \cdot (${est0_str}) + w_2 \cdot (${G_s1_str}) \cdot (${est1_str}) \\
            &= (${w0_str}) \cdot (${G_est0_str}) + (${w1_str}) \cdot (${G_est1_str}) \\
            &= (${steady_weighted0_str}) + (${steady_weighted1_str}) = ${y_steady_real.toFixed(4)} \quad \text{bei} \quad t = ${t.toFixed(2)}
            \end{aligned}`,
            document.getElementById('steady-state-response-now-eq'),
            { displayMode: true, throwOnError: false }
        );
        
    } catch (error) {
        console.error('Error rendering time-dependent equations:', error);
    }
}

// Render all equations (static and time-dependent)
function renderEquations() {
    renderEquationsStatic();
    renderEquationsTimeDependent();
}

function resimulate() {
    // Update state-space matrices from current parameters
    state.matrices = system.getStateSpaceMatrices();
    
    // Compute complex frequency parameters for sinusoidal input u(t) = A*sin(2πft)
    // Using Euler's formula: sin(ωt) = (1/2j)[e^(jωt) - e^(-jωt)]
    // Therefore: A*sin(ωt) = A*(-j/2)*e^(jωt) + A*(j/2)*e^(-jωt)
    const omega = 2 * Math.PI * state.f;
    const A = state.A;
    state.s_values = [
        { s: math.complex(0, omega),   weight: math.complex(0, -0.5 * A) },  // s = +jω, weight = -jA/2
        { s: math.complex(0, -omega),  weight: math.complex(0, 0.5 * A) }    // s = -jω, weight = +jA/2
    ];
    
    // Compute transfer function data (parameter-dependent but time-independent)
    // This makes the data available for both simulation and equation rendering
    system.computeTransferFunctionData();
    
    // Recompute Bode diagram with updated parameters
    if (bodeRenderer) {
        bodeRenderer.computeFrequencyResponse();
    }
    
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
    document.getElementById('m0-slider').addEventListener('input', (e) => {
        state.m0 = parseFloat(e.target.value);
        document.getElementById('m0-value').textContent = state.m0.toFixed(2);
        resimulate();
    });
    
    document.getElementById('m1-slider').addEventListener('input', (e) => {
        state.m1 = parseFloat(e.target.value);
        document.getElementById('m1-value').textContent = state.m1.toFixed(2);
        resimulate();
    });
    
    document.getElementById('c0-slider').addEventListener('input', (e) => {
        state.c0 = parseFloat(e.target.value);
        document.getElementById('c0-value').textContent = state.c0.toFixed(1);
        resimulate();
    });
    
    document.getElementById('c1-slider').addEventListener('input', (e) => {
        state.c1 = parseFloat(e.target.value);
        document.getElementById('c1-value').textContent = state.c1.toFixed(1);
        resimulate();
    });
    
    document.getElementById('c2-slider').addEventListener('input', (e) => {
        state.c2 = parseFloat(e.target.value);
        document.getElementById('c2-value').textContent = state.c2.toFixed(1);
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
        renderEquationsTimeDependent();  // Update time-dependent equations
        render();
    });
    
    // Offset slider
    document.getElementById('offset-slider').addEventListener('input', (e) => {
        state.offset = parseFloat(e.target.value);
        document.getElementById('offset-value').textContent = state.offset.toFixed(0);
        render();
    });
    
    // Trajectory display checkboxes
    document.getElementById('show-input').addEventListener('change', (e) => {
        state.isShowingInput = e.target.checked;
        render();
    });
    
    document.getElementById('show-analytic-total').addEventListener('change', (e) => {
        state.isShowingAnalyticTotal = e.target.checked;
        render();
    });
    
    document.getElementById('show-analytic-transient').addEventListener('change', (e) => {
        state.isShowingAnalyticTransient = e.target.checked;
        render();
    });
    
    document.getElementById('show-analytic-steady').addEventListener('change', (e) => {
        state.isShowingAnalyticSteady = e.target.checked;
        render();
    });
    
    document.getElementById('show-numerical').addEventListener('change', (e) => {
        state.isShowingNumerical = e.target.checked;
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
            bodeRenderer.setupCanvas();
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
    
    // Enable auto-rendering of inline math with \(...\) delimiters
    renderMathInElement(document.body, {
        delimiters: [
            {left: "\\(", right: "\\)", display: false}
        ],
        throwOnError: false
    });
});
