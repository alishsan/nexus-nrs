// DWBA Web Dashboard JavaScript

/** Default channel partial waves for Numerov/phase, total σ (cross-section tab), inelastic, elastic metadata. */
const DEFAULT_CHANNEL_L = '0,1,2,3,4,5';

/** Default lab energies (MeV) for Phase Shifts and R-Matrices tabs only. */
const DEFAULT_PHASE_RMATRIX_ENERGIES = '4,8,12,16,20';

/** Plotly: smooth curves between samples (cubic spline). `smoothing` in [0, 1.3]; moderate value limits overshoot on oscillatory DCS. */
const PLOTLY_SPLINE_LINE = { shape: 'spline', smoothing: 0.65 };

/** Sort scatter points by x so spline and line traces are well-behaved. */
function sortTraceByX(trace) {
    if (!trace || !trace.x || !trace.y || trace.x.length !== trace.y.length) return;
    const pairs = trace.x.map((x, i) => [Number(x), trace.y[i]]);
    pairs.sort((a, b) => a[0] - b[0]);
    trace.x = pairs.map(p => p[0]);
    trace.y = pairs.map(p => p[1]);
}

class DWBADashboard {
    constructor() {
        // Use same-origin for API when served over http/https (e.g. Railway); use localhost only when opened from file://
        const dashboardUrl = 'http://localhost:3000';
        const apiFromQuery = new URLSearchParams(window.location.search).get('api');
        const sameOrigin = window.location.protocol === 'http:' || window.location.protocol === 'https:';
        this.apiBase = apiFromQuery !== null ? apiFromQuery : (sameOrigin ? '' : dashboardUrl);
        this.currentData = null;
        this.initializeEventListeners();
        this.loadDefaultParameters();
        this.checkApiHealth();
        this.loadTransferDefault();  // default DCS on Transfer tab (16O d-p preset when selected)
        if (this.apiBase === dashboardUrl) {
            this.showApiNotice(dashboardUrl);
        }
    }

    /** Load default DCS for the Transfer tab (`/api/transfer-default`). */
    async loadTransferDefault() {
        const targetEl = document.getElementById('transfer_target');
        const target = targetEl ? targetEl.value : '16O';
        const rtEl = document.getElementById('reaction_type');
        const reactionType = rtEl ? rtEl.value : 'p-d';
        const labelEl = document.getElementById('transfer-default-label');
        const q = new URLSearchParams({ target, reaction_type: reactionType });
        if (target === 'generic') {
            const aEl = document.getElementById('transfer_target_A');
            const zEl = document.getElementById('transfer_target_Z');
            if (aEl && aEl.value !== '') q.set('transfer_target_A', String(parseInt(aEl.value, 10) || 40));
            if (zEl && zEl.value !== '') q.set('transfer_target_Z', String(parseInt(zEl.value, 10) || 20));
        }
        try {
            const r = await fetch(`${this.apiBase}/api/transfer-default?${q.toString()}`);
            if (!r.ok) {
                if (labelEl) labelEl.textContent = `${target} ${reactionType} (load failed)`;
                return;
            }
            const json = await r.json();
            if (json.success && json.data && json.data.transfer && json.data.transfer.length) {
                this.currentData = this.currentData || {};
                this.currentData.transfer = json.data.transfer;
                this.currentData.transferTargetLabel = json.target || json.data?.parameters?.target || target;
                if (labelEl) labelEl.textContent = this.currentData.transferTargetLabel;
                this.plotTransfer();
            } else if (json.error && labelEl) {
                labelEl.textContent = `${target} ${reactionType} (error)`;
            }
        } catch (e) {
            console.warn('Transfer default not loaded:', e.message);
            if (labelEl) labelEl.textContent = `${target} ${reactionType} (offline)`;
        }
    }

    showApiNotice(dashboardUrl) {
        const statusDiv = document.getElementById('status-messages');
        if (!statusDiv) return;
        statusDiv.innerHTML = '<div class="alert alert-info mb-0" role="alert"><strong>Calculations use the dashboard server.</strong> Start it with: <code>cd web-dashboard && lein run</code>. Then <a href="' + dashboardUrl + '" class="alert-link">open this app from ' + dashboardUrl + '</a> so the page and API come from the same server.</div>';
    }

    async checkApiHealth() {
        try {
            const r = await fetch(`${this.apiBase}/api/health`);
            if (!r.ok) throw new Error('API not OK');
        } catch (e) {
            const statusDiv = document.getElementById('status-messages');
            if (statusDiv) {
                const isLocalDev = this.apiBase === 'http://localhost:3000';
                const msg = isLocalDev
                    ? '<div class="error"><i class="fas fa-exclamation-triangle"></i> Dashboard API not reachable. Start the server: <code>cd web-dashboard && lein run</code>, then <a href="http://localhost:3000">open http://localhost:3000</a> in your browser.</div>'
                    : '<div class="error"><i class="fas fa-exclamation-triangle"></i> Dashboard API not reachable. Please try again later or check the server.</div>';
                statusDiv.innerHTML = msg;
            }
        }
    }

    initializeEventListeners() {
        // Parameter sliders
        ['V0', 'R0', 'a0', 'radius'].forEach(param => {
            const slider = document.getElementById(param);
            const valueDisplay = document.getElementById(`${param}-value`);
            
            slider.addEventListener('input', (e) => {
                const value = parseFloat(e.target.value);
                const unit = param === 'V0' ? 'MeV' : 'fm';
                valueDisplay.textContent = `${value} ${unit}`;
            });
        });

        // Reset button
        document.getElementById('reset-btn')?.addEventListener('click', () => {
            this.resetParameters();
        });

        // Calculate buttons (one per tab) – bound in JS so no globals needed
        const calcButtons = [
            ['calculate-phase-btn', 'calculatePhase'],
            ['calculate-rmatrix-btn', 'calculateRMatrix'],
            ['calculate-potential-btn', 'calculatePotentials'],
            ['calculate-cross-section-btn', 'calculateCrossSections'],
            ['calculate-elastic-btn', 'calculateElastic'],
            ['calculate-inelastic-btn', 'calculateInelastic'],
            ['calculate-transfer-btn', 'calculateTransfer'],
            ['calculate-dashboard-btn', 'calculateDashboard']
        ];
        calcButtons.forEach(([id, method]) => {
            const btn = document.getElementById(id);
            if (btn) btn.addEventListener('click', () => this[method]());
        });

        // Target / reaction: reload default plot
        document.getElementById('transfer_target')?.addEventListener('change', () => {
            this.toggleTransferGenericInputs();
            this.loadTransferDefault();
        });
        this.toggleTransferGenericInputs();
        ['transfer_target_A', 'transfer_target_Z'].forEach(id => {
            document.getElementById(id)?.addEventListener('change', () => {
                if ((document.getElementById('transfer_target') || {}).value === 'generic') {
                    this.loadTransferDefault();
                }
            });
        });
        document.getElementById('reaction_type')?.addEventListener('change', () => {
            this.loadTransferDefault();
        });

        // Inelastic target: set E_ex and β from preset
        const inelasticPresets = {
            '12C': { E_ex: 4.44, beta: 0.25, lambda: '2' },
            '16O': { E_ex: 6.13, beta: 0.2, lambda: '2' }
        };
        document.getElementById('inelastic_target')?.addEventListener('change', () => {
            const sel = document.getElementById('inelastic_target');
            const preset = sel && inelasticPresets[sel.value];
            if (preset) {
                const eEx = document.getElementById('E_ex');
                const beta = document.getElementById('beta');
                const lambda = document.getElementById('lambda');
                if (eEx) eEx.value = preset.E_ex;
                if (beta) beta.value = preset.beta;
                if (lambda) lambda.value = preset.lambda;
            }
        });

        const syncInelasticModelUi = () => {
            const m = (document.getElementById('inelastic_model') || {}).value || 'standard';
            const li11Box = document.getElementById('li11_paper_options');
            const paperish = m === 'li11_paper1708' || m === 'tanaka2017';
            if (li11Box) li11Box.style.display = paperish ? '' : 'none';
        };
        document.getElementById('inelastic_model')?.addEventListener('change', () => {
            syncInelasticModelUi();
            const m = (document.getElementById('inelastic_model') || {}).value || 'standard';
            if (m === 'li11_paper1708' || m === 'tanaka2017') {
                const eEx = document.getElementById('E_ex');
                if (eEx && (parseFloat(eEx.value) === 4.44 || parseFloat(eEx.value) === 6.13)) eEx.value = '0.8';
                const proj = document.getElementById('inelastic_projectile');
                if (proj) proj.value = 'p';
            }
        });
        syncInelasticModelUi();

        // Elastic target: show A/Z inputs only when Generic is selected
        document.getElementById('elastic_target')?.addEventListener('change', () => this.toggleElasticGenericInputs());
        this.toggleElasticGenericInputs();

        document.getElementById('elastic_ratio_rutherford')?.addEventListener('change', () => {
            if (this.currentData?.elastic) this.plotElastic();
        });
    }

    /** Lab projectile mass (MeV/c²) and charge number for elastic kinematics (aligned with backend). */
    elasticProjectileKinematics() {
        const proj = (document.getElementById('elastic_projectile') || {}).value || 'p';
        const masses = { p: 938.272, n: 939.565, d: 1875.613, a: 3727.379 };
        const Z = { p: 1, n: 0, d: 1, a: 2 };
        return { mass: masses[proj] || 938.272, Z: Z[proj] !== undefined ? Z[proj] : 1 };
    }

    /** Target A, Z for elastic tab (aligned with backend presets). */
    elasticTargetAZ() {
        const tgt = (document.getElementById('elastic_target') || {}).value || '16O';
        if (tgt === '12C') return { A: 12, Z: 6 };
        if (tgt === '16O') return { A: 16, Z: 8 };
        if (tgt === '148Sm') return { A: 148, Z: 62 };
        const A = parseInt((document.getElementById('elastic_target_A') || {}).value, 10);
        const Z = parseInt((document.getElementById('elastic_target_Z') || {}).value, 10);
        return {
            A: (Number.isFinite(A) && A >= 1) ? A : 16,
            Z: (Number.isFinite(Z) && Z >= 0) ? Z : 8
        };
    }

    /**
     * Rutherford dσ/dΩ in mb/sr (non-relativistic, point Coulomb, CM frame).
     * E_cm_mev = lab projectile kinetic × m_target / (m_proj + m_target); θ_cm in degrees.
     */
    rutherfordDsDomegaMbSr(Z1, Z2, E_cm_MeV, thetaCmDeg) {
        const e2 = 1.44; // MeV·fm
        if (Z1 * Z2 <= 0 || !(E_cm_MeV > 0)) return null;
        const th = (thetaCmDeg * Math.PI) / 180;
        const s = Math.sin(th / 2);
        if (s < 1e-5) return null;
        // dσ/dΩ (mb/sr): Coulomb prefactor / sin⁴(θ/2), conventional ×10 from fm²-area units
        return Math.pow((Z1 * Z2 * e2) / (4 * E_cm_MeV), 2) / Math.pow(s, 4) * 10.0;
    }

    labToCmKineticMeV(mProjMeV, mTargMeV, eLabMeV) {
        return eLabMeV * (mTargMeV / (mProjMeV + mTargMeV));
    }

    toggleElasticGenericInputs() {
        const targetEl = document.getElementById('elastic_target');
        const show = targetEl && targetEl.value === 'generic';
        const rowA = document.getElementById('elastic-generic-AZ-row');
        const rowZ = document.getElementById('elastic-generic-Z-row');
        if (rowA) rowA.style.display = show ? '' : 'none';
        if (rowZ) rowZ.style.display = show ? '' : 'none';
    }

    toggleTransferGenericInputs() {
        const targetEl = document.getElementById('transfer_target');
        const show = targetEl && targetEl.value === 'generic';
        const row = document.getElementById('transfer-generic-AZ-row');
        if (row) row.style.display = show ? '' : 'none';
    }

    _setButtonLoading(btnId, loading) {
        const btn = document.getElementById(btnId);
        if (!btn) return;
        btn.disabled = loading;
        btn.innerHTML = loading ? '<i class="fas fa-spinner fa-spin"></i> Calculating...' : '<i class="fas fa-calculator"></i> Calculate';
    }

    async _post(url, params) {
        const fullUrl = `${this.apiBase}${url}`;
        const response = await fetch(fullUrl, {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify(params)
        });
        const text = await response.text();
        if (!response.ok) {
            throw new Error(`HTTP ${response.status}: ${text.slice(0, 200)}`);
        }
        let result;
        try {
            result = JSON.parse(text);
        } catch (e) {
            throw new Error('Invalid JSON response: ' + text.slice(0, 200));
        }
        return result;
    }

    /** Normalize API response so we always have underscore keys (phase_shifts, r_matrices, etc.) */
    _normalizeData(data) {
        if (!data) return data;
        return {
            phase_shifts: data.phase_shifts || data.phaseShifts || data['phase-shifts'],
            r_matrices: data.r_matrices || data.rMatrices || data['r-matrices'],
            potentials: data.potentials,
            cross_sections: data.cross_sections || data.crossSections || data['cross-sections'],
            parameters: data.parameters,
            elastic: data.elastic,
            inelastic: data.inelastic,
            transfer: data.transfer
        };
    }

    async loadDefaultParameters() {
        try {
            const response = await fetch(`${this.apiBase}/api/parameters`);
            const data = await response.json();
            
            if (data.default_parameters) {
                this.setParameters(data.default_parameters);
            }
        } catch (error) {
            console.error('Error loading default parameters:', error);
        }
    }

    setParameters(params) {
        const setIf = (id, v) => {
            const el = document.getElementById(id);
            if (el && v !== undefined && v !== null) el.value = v;
        };
        document.getElementById('V0').value = params.V0;
        document.getElementById('R0').value = params.R0;
        document.getElementById('a0').value = params.a0;
        document.getElementById('radius').value = params.radius;
        const eVal = Array.isArray(params.energies) ? params.energies.join(',') : (params.energies ?? '50');
        const lVal = Array.isArray(params.L_values) ? params.L_values.join(',') : (params.L_values ?? DEFAULT_CHANNEL_L);
        const setVal = (id, v) => { const el = document.getElementById(id); if (el) el.value = v; };
        const phaseE = Array.isArray(params.phase_energies) ? params.phase_energies.join(',') : (params.phase_energies != null ? String(params.phase_energies) : DEFAULT_PHASE_RMATRIX_ENERGIES);
        const rmatE = Array.isArray(params.rmatrix_energies) ? params.rmatrix_energies.join(',') : (params.rmatrix_energies != null ? String(params.rmatrix_energies) : DEFAULT_PHASE_RMATRIX_ENERGIES);
        setVal('phase-energy-range', phaseE);
        setVal('rmatrix-energy-range', rmatE);
        ['potential-energy-range', 'dashboard-energy-range', 'cross-section-energy-range', 'elastic-lab-energies', 'inelastic-beam-energies'].forEach(id => setVal(id, eVal));
        ['phase-L-values', 'rmatrix-L-values', 'potential-L-values', 'dashboard-L-values'].forEach(id => setVal(id, lVal));
        
        if (params.E_ex !== undefined) document.getElementById('E_ex').value = params.E_ex;
        if (params.lambdas !== undefined) document.getElementById('lambda').value = Array.isArray(params.lambdas) ? params.lambdas.join(',') : String(params.lambdas);
        else if (params.lambda !== undefined) document.getElementById('lambda').value = params.lambda;
        if (params.beta !== undefined) document.getElementById('beta').value = params.beta;
        if (params.reaction_type !== undefined) document.getElementById('reaction_type').value = params.reaction_type;
        if (params.elastic_projectile !== undefined && document.getElementById('elastic_projectile')) {
            document.getElementById('elastic_projectile').value = params.elastic_projectile;
        }
        if (params.elastic_target !== undefined && document.getElementById('elastic_target')) {
            document.getElementById('elastic_target').value = params.elastic_target;
        }
        if (params.elastic_target_A !== undefined) {
            const el = document.getElementById('elastic_target_A');
            if (el) el.value = params.elastic_target_A;
        }
        if (params.elastic_target_Z !== undefined) {
            const el = document.getElementById('elastic_target_Z');
            if (el) el.value = params.elastic_target_Z;
        }
        if (params.elastic_dsigma_L_max !== undefined) {
            const el = document.getElementById('elastic_dsigma_L_max');
            if (el) el.value = params.elastic_dsigma_L_max;
        }
        setIf('elastic_W0', params.W0);
        setIf('elastic_RW', params.R_W);
        setIf('elastic_aW', params.a_W);
        this.toggleElasticGenericInputs && this.toggleElasticGenericInputs();

        // Update slider displays
        document.getElementById('V0-value').textContent = `${params.V0} MeV`;
        document.getElementById('R0-value').textContent = `${params.R0} fm`;
        document.getElementById('a0-value').textContent = `${params.a0} fm`;
        document.getElementById('radius-value').textContent = `${params.radius} fm`;

        // Transfer-tab defaults (optional keys from /api/parameters)
        setIf('transfer_V0', params.V0);
        setIf('transfer_R0', params.R0);
        setIf('transfer_a0', params.a0);
        setIf('transfer_W0', params.transfer_W0);
        setIf('transfer_RW', params.transfer_RW ?? params.R_W);
        setIf('transfer_aW', params.transfer_aW ?? params.a_W);
        setIf('transfer_r_max', params.transfer_r_max);
        setIf('transfer_h', params.transfer_h);
        setIf('transfer_S', params.transfer_S);
        setIf('transfer_partial_L', params.transfer_partial_L);
        setIf('transfer_l_i', params.transfer_l_i);
        setIf('transfer_l_f', params.transfer_l_f);
        if (params.transfer_energies !== undefined) {
            const te = document.getElementById('transfer_energies');
            if (te) te.value = Array.isArray(params.transfer_energies) ? params.transfer_energies.join(',') : String(params.transfer_energies);
        }
        if (params.transfer_L_values !== undefined) {
            const tl = document.getElementById('transfer_L_values');
            if (tl) tl.value = Array.isArray(params.transfer_L_values) ? params.transfer_L_values.join(',') : String(params.transfer_L_values);
        }
    }

    /**
     * Which **E** and **L** inputs to use for `/api/calculate` (phase, R-matrix, potential, total σ, dashboard).
     * Cross-section tab has **energy only** in the UI; **L** defaults to {@link DEFAULT_CHANNEL_L} for the server sum.
     * @param {'phase'|'rmatrix'|'potential'|'cross'|'dashboard'} tabKey
     */
    getCoreEnergiesL(tabKey) {
        const defaultE = '50';
        const pairs = {
            phase: ['phase-energy-range', 'phase-L-values'],
            rmatrix: ['rmatrix-energy-range', 'rmatrix-L-values'],
            potential: ['potential-energy-range', 'potential-L-values'],
            cross: ['cross-section-energy-range', null],
            dashboard: ['dashboard-energy-range', 'dashboard-L-values'],
        };
        const [eId, lId] = pairs[tabKey] || pairs.phase;
        const rawE = ((document.getElementById(eId) || {}).value || '').trim();
        const fallbackE = (tabKey === 'phase' || tabKey === 'rmatrix') ? DEFAULT_PHASE_RMATRIX_ENERGIES : defaultE;
        const eStr = rawE || fallbackE;
        const lStr = lId
            ? (((document.getElementById(lId) || {}).value || '').trim() || DEFAULT_CHANNEL_L)
            : DEFAULT_CHANNEL_L;
        return { energies: eStr, L_values: lStr };
    }

    /** Map Calculate button id → tab key for {@link getCoreEnergiesL}. */
    coreTabKeyFromButton(btnId) {
        const m = {
            'calculate-phase-btn': 'phase',
            'calculate-rmatrix-btn': 'rmatrix',
            'calculate-potential-btn': 'potential',
            'calculate-cross-section-btn': 'cross',
            'calculate-dashboard-btn': 'dashboard',
        };
        return m[btnId] || 'phase';
    }

    getParameters() {
        const num = (id, fallback) => {
            const el = document.getElementById(id);
            const v = el ? parseFloat(el.value) : NaN;
            return (v !== undefined && !Number.isNaN(v)) ? v : fallback;
        };
        const int = (id, fallback) => {
            const el = document.getElementById(id);
            const v = el ? parseInt(el.value, 10) : NaN;
            return (v !== undefined && !Number.isNaN(v)) ? v : fallback;
        };
        // Default **phase** tab fields — used for generic API payloads; per-tab `getCoreEnergiesL` overrides for /api/calculate
        const phaseE = (document.getElementById('phase-energy-range') || {}).value || DEFAULT_PHASE_RMATRIX_ENERGIES;
        const phaseL = (document.getElementById('phase-L-values') || {}).value || DEFAULT_CHANNEL_L;
        return {
            V0: num('V0', 65),
            R0: num('R0', 7.5),
            a0: num('a0', 0.67),
            radius: num('radius', 100),
            energies: String(phaseE).trim() || DEFAULT_PHASE_RMATRIX_ENERGIES,
            L_values: String(phaseL).trim() || DEFAULT_CHANNEL_L,
            E_ex: num('E_ex', 4.44),
            lambda: int('lambda', 2),
            lambdas: (document.getElementById('lambda') || {}).value || '2',  // e.g. "2" or "2,3,4" for multiple multipoles
            beta: num('beta', 0.25),
            reaction_type: (document.getElementById('reaction_type') || {}).value || 'p-d',
            target: (document.getElementById('transfer_target') || {}).value || '16O',
            projectile: (document.getElementById('inelastic_projectile') || {}).value || 'p',
            inelastic_target: (document.getElementById('inelastic_target') || {}).value || '12C',
            elastic_projectile: (document.getElementById('elastic_projectile') || {}).value || 'a',
            elastic_target: (document.getElementById('elastic_target') || {}).value || '148Sm',
            elastic_target_A: int('elastic_target_A', 148),
            elastic_target_Z: int('elastic_target_Z', 62),
            elastic_dsigma_L_max: int('elastic_dsigma_L_max', 41),
            // Complex Woods-Saxon for elastic → **`functions/*elastic-imag-ws-params*`**
            W0: num('elastic_W0', 30),
            R_W: num('elastic_RW', 7.5),
            a_W: num('elastic_aW', 0.67),
            // Transfer tab (POST /api/transfer)
            transfer_energies: ((document.getElementById('transfer_energies') || {}).value || '').trim(),
            transfer_L_values: (document.getElementById('transfer_L_values') || {}).value || '1',
            transfer_S: num('transfer_S', 1),
            transfer_partial_L: int('transfer_partial_L', 1),
            transfer_l_i: int('transfer_l_i', 1),
            transfer_l_f: int('transfer_l_f', 0),
            transfer_V0: num('transfer_V0', 40),
            transfer_R0: num('transfer_R0', 2),
            transfer_a0: num('transfer_a0', 0.6),
            transfer_W0: num('transfer_W0', 0),
            transfer_RW: num('transfer_RW', 2),
            transfer_aW: num('transfer_aW', 0.6),
            transfer_r_max: num('transfer_r_max', 20),
            transfer_h: num('transfer_h', 0.01),
            transfer_target_A: int('transfer_target_A', 40),
            transfer_target_Z: int('transfer_target_Z', 20),
            inelastic_model: (document.getElementById('inelastic_model') || {}).value || 'standard',
            li11_optical_set: (document.getElementById('li11_optical_set') || {}).value || 'V',
            li11_beta_scale: num('li11_beta_scale', 1),
            li11_L_max: int('li11_L_max', 22),
            li11_r_max: num('li11_r_max', 45),
            li11_h: num('li11_h', 0.02),
            li11_theta_step: num('li11_theta_step', 5),
            li11_transfer_ell: int('li11_transfer_ell', 1)
        };
    }

    showStatus(message, type = 'info') {
        const statusDiv = document.getElementById('status-messages');
        const alertClass = type === 'error' ? 'error' : type === 'success' ? 'success' : 'alert-info';
        
        statusDiv.innerHTML = `
            <div class="${alertClass}">
                <i class="fas fa-${type === 'error' ? 'exclamation-triangle' : type === 'success' ? 'check-circle' : 'info-circle'}"></i>
                ${message}
            </div>
        `;
        
        // Auto-hide after 8 seconds (longer so user sees success/error)
        setTimeout(() => {
            statusDiv.innerHTML = '';
        }, 8000);
    }

    async _runCoreCalculation(btnId) {
        const tabKey = this.coreTabKeyFromButton(btnId);
        const { energies, L_values } = this.getCoreEnergiesL(tabKey);
        const energiesStr = String(energies ?? '').trim();
        const LValuesStr = String(L_values ?? '').trim();
        if (!energiesStr) {
            this.showStatus('Please set an energy range (MeV) on this tab.', 'error');
            return;
        }
        if (tabKey !== 'cross' && !LValuesStr) {
            this.showStatus('Please set channel partial waves L on this tab.', 'error');
            return;
        }
        const params = { ...this.getParameters(), energies: energiesStr, L_values: LValuesStr };
        this._setButtonLoading(btnId, true);
        this.showStatus('Calculating...', 'info');
        const startTime = Date.now();
        try {
            const result = await this._post('/api/calculate', params);
            if (!result.success) throw new Error(result.error || 'Calculation failed');
            this.currentData = this._normalizeData(result.data);
            if (!this.currentData) throw new Error('No data in response');
            const n = (this.currentData.phase_shifts || []).length;
            console.log('Core calculation OK, data keys:', Object.keys(this.currentData), 'phase_shifts count:', n);
            if (n === 0) {
                this.showStatus('Calculation returned no data points. Check energy range and L on this tab.', 'error');
                return;
            }
            try {
                this.updateAllPlots();
                this.updateDashboardStats(Date.now() - startTime);
                this.showStatus(`Done in ${Date.now() - startTime}ms — ${n} points plotted.`, 'success');
            } catch (plotError) {
                console.error('Plot error:', plotError);
                this.showStatus(`Data received but plot failed: ${plotError.message}`, 'error');
            }
        } catch (error) {
            console.error('Calculation error:', error);
            let msg = error.message;
            if (msg.includes('404')) {
                msg = this.apiBase === 'http://localhost:3000'
                    ? 'API not found. Start the server: cd web-dashboard && lein run — then open http://localhost:3000 in your browser (click the link above).'
                    : 'API not found. The server may be starting or temporarily unavailable.';
            }
            this.showStatus(`Error: ${msg}`, 'error');
        } finally {
            this._setButtonLoading(btnId, false);
        }
    }

    async _runTabCalculation(btnId, path, dataKey) {
        const params = this.getParameters();
        const energiesStr = String(params.energies ?? '').trim();
        const LValuesStr = String(params.L_values ?? '').trim();
        if (!energiesStr || !LValuesStr) {
            this.showStatus('Please provide valid energy range and angular momenta', 'error');
            return;
        }
        this._setButtonLoading(btnId, true);
        this.showStatus('Calculating...', 'info');
        const startTime = Date.now();
        try {
            const result = await this._post(path, params);
            if (!result.success) throw new Error(result.error || 'Calculation failed');
            this.currentData = this.currentData || {};
            const raw = result.data && (result.data[dataKey] || result.data[dataKey.replace(/_/g, '-')]);
            if (raw) this.currentData[dataKey] = raw;
            this.updateAllPlots();
            if (dataKey === 'elastic' || dataKey === 'inelastic' || dataKey === 'transfer') {
                this.showStatus(`Done in ${Date.now() - startTime}ms`, 'success');
            } else {
                this.updateDashboardStats(Date.now() - startTime);
                this.showStatus(`Done in ${Date.now() - startTime}ms`, 'success');
            }
        } catch (error) {
            console.error('Calculation error:', error);
            let msg = error.message;
            if (msg.includes('404')) {
                msg = this.apiBase === 'http://localhost:3000'
                    ? 'API not found. Start the server: cd web-dashboard && lein run — then open http://localhost:3000 in your browser (click the link above).'
                    : 'API not found. The server may be starting or temporarily unavailable.';
            }
            this.showStatus(`Error: ${msg}`, 'error');
        } finally {
            this._setButtonLoading(btnId, false);
        }
    }

    async calculatePhase() { await this._runCoreCalculation('calculate-phase-btn'); }
    async calculateRMatrix() { await this._runCoreCalculation('calculate-rmatrix-btn'); }
    async calculatePotentials() { await this._runCoreCalculation('calculate-potential-btn'); }
    async calculateCrossSections() { await this._runCoreCalculation('calculate-cross-section-btn'); }
    async calculateDashboard() { await this._runCoreCalculation('calculate-dashboard-btn'); }
    async calculateElastic() {
        const energyEl = document.getElementById('elastic-lab-energies');
        const energiesStr = (energyEl && energyEl.value || '').trim();
        if (!energiesStr) {
            this.showStatus('Please set lab energies (MeV) on the Elastic tab.', 'error');
            return;
        }
        const params = this.getParameters();
        const body = { ...params, energies: energiesStr, L_values: DEFAULT_CHANNEL_L };
        
        // Debug: log what we're sending
        console.log('Elastic calculation - sending energies:', energiesStr);
        
        this._setButtonLoading('calculate-elastic-btn', true);
        this.showStatus('Calculating...', 'info');
        const startTime = Date.now();
        try {
            const result = await this._post('/api/elastic', body);
            if (!result.success) throw new Error(result.error || 'Calculation failed');
            
            // Debug: log what we received
            const receivedEnergies = result.data && result.data.parameters && result.data.parameters.energies;
            console.log('Elastic calculation - received energies:', receivedEnergies);
            
            this.currentData = this.currentData || {};
            const raw = result.data && (result.data.elastic || result.data['elastic']);
            // Replace elastic data completely (don't merge with old data)
            if (raw) {
                this.currentData.elastic = raw;
                // Debug: log unique energies in the data
                const uniqueEnergies = [...new Set(raw.map(p => p.energy))].sort((a, b) => a - b);
                console.log('Elastic calculation - unique energies in plot data:', uniqueEnergies);
            } else {
                delete this.currentData.elastic;
            }
            this.updateAllPlots();
            const used = result.data && result.data.parameters && result.data.parameters.energies;
            const n = Array.isArray(used) ? used.length : 0;
            const p = result.data && result.data.parameters;
            const edL = p && p.elastic_dsigma_L_values;
            let lPart = '';
            if (Array.isArray(edL) && edL.length) {
                const lo = edL[0];
                const hi = edL[edL.length - 1];
                lPart = `, elastic dσ: L=${lo}…${hi} (${edL.length} partial waves)`;
            } else if (p && Array.isArray(p.L_values) && p.L_values.length) {
                lPart = `, L=[${p.L_values.join(',')}]`;
            }
            this.showStatus(`Done in ${Date.now() - startTime}ms${n ? ` (${n} energies)` : ''}${lPart}`, 'success');
        } catch (error) {
            console.error('Elastic calculation error:', error);
            this.showStatus(`Error: ${error.message}`, 'error');
        } finally {
            this._setButtonLoading('calculate-elastic-btn', false);
        }
    }
    async calculateInelastic() {
        const modelRaw = String((document.getElementById('inelastic_model') || {}).value || 'standard').trim();
        const model = modelRaw.toLowerCase().replace(/-/g, '_');
        const li11Models = new Set(['li11_paper1708', 'li11', 'tanaka2017', 'paper1708']);
        const isLi11Paper = li11Models.has(model);

        const energyInel = document.getElementById('inelastic-beam-energies');
        const energiesStr = (energyInel && energyInel.value || '').trim();

        if (!energiesStr) {
            this.showStatus('Please set beam energies (MeV) on the Inelastic tab. For ¹¹Li paper mode this is E_p,lab.', 'error');
            return;
        }

        const params = { ...this.getParameters(), energies: energiesStr };
        params.L_values = isLi11Paper ? '0' : DEFAULT_CHANNEL_L;

        this._setButtonLoading('calculate-inelastic-btn', true);
        this.showStatus('Calculating...', 'info');
        const startTime = Date.now();
        try {
            const result = await this._post('/api/inelastic', params);
            if (!result.success) throw new Error(result.error || 'Calculation failed');
            this.currentData = this.currentData || {};
            const raw = result.data && (result.data.inelastic || result.data['inelastic']);
            if (raw) this.currentData.inelastic = raw;
            this.updateAllPlots();
            this.showStatus(`Done in ${Date.now() - startTime}ms`, 'success');
        } catch (error) {
            console.error('Calculation error:', error);
            let msg = error.message;
            if (msg.includes('404')) {
                msg = this.apiBase === 'http://localhost:3000'
                    ? 'API not found. Start the server: cd web-dashboard && lein run — then open http://localhost:3000 in your browser (click the link above).'
                    : 'API not found. The server may be starting or temporarily unavailable.';
            }
            this.showStatus(`Error: ${msg}`, 'error');
        } finally {
            this._setButtonLoading('calculate-inelastic-btn', false);
        }
    }

    /** Transfer: prefer energies / L from this tab; fall back to main panel; default energy 20 MeV if both empty. */
    async calculateTransfer() {
        const tE = (document.getElementById('transfer_energies')?.value || '').trim();
        const mainE = (document.getElementById('phase-energy-range')?.value || '').trim();
        const energiesStr = tE || mainE || '20';
        const tL = (document.getElementById('transfer_L_values')?.value || '').trim();
        const LValuesStr = tL || '1';
        const params = { ...this.getParameters(), energies: energiesStr, L_values: LValuesStr };
        if (tE) params.transfer_energies = tE;

        this._setButtonLoading('calculate-transfer-btn', true);
        this.showStatus('Calculating transfer...', 'info');
        const startTime = Date.now();
        try {
            const result = await this._post('/api/transfer', params);
            if (!result.success) throw new Error(result.error || 'Calculation failed');
            this.currentData = this.currentData || {};
            const raw = result.data && (result.data.transfer || result.data['transfer']);
            if (raw) this.currentData.transfer = raw;
            const lab = result.data?.parameters?.reaction_type;
            const s = result.data?.parameters?.S_factor;
            const impl = result.data?.parameters?.implementation;
            const lbl = document.getElementById('transfer-default-label');
            if (lbl && lab != null && s != null) {
                lbl.textContent = impl ? `${lab}, S=${s} — ${impl}` : `${lab}, S=${s} (calculated)`;
            }
            this.updateAllPlots();
            const implShort = impl && impl.length > 42 ? impl.substring(0, 40) + '…' : impl;
            this.showStatus(`Transfer done in ${Date.now() - startTime}ms${implShort ? ` (${implShort})` : ''}`, 'success');
        } catch (error) {
            console.error('Calculation error:', error);
            let msg = error.message;
            if (msg.includes('404')) {
                msg = this.apiBase === 'http://localhost:3000'
                    ? 'API not found. Start the server: cd web-dashboard && lein run — then open http://localhost:3000 in your browser (click the link above).'
                    : 'API not found. The server may be starting or temporarily unavailable.';
            }
            this.showStatus(`Error: ${msg}`, 'error');
        } finally {
            this._setButtonLoading('calculate-transfer-btn', false);
        }
    }

    updateAllPlots() {
        if (!this.currentData) return;

        this.plotPhaseShifts();
        this.plotRMatrices();
        this.plotPotentials();
        this.plotCrossSections();
        if (this.currentData.elastic) this.plotElastic();
        if (this.currentData.inelastic) this.plotInelastic();
        if (this.currentData.transfer) this.plotTransfer();
        this.plotDashboard();
    }

    plotPhaseShifts() {
        const data = this.currentData.phase_shifts || this.currentData['phase-shifts'] || this.currentData.phaseShifts;
        if (!data || !Array.isArray(data) || data.length === 0) {
            console.warn('plotPhaseShifts: no data', { hasData: !!data, length: data?.length });
            return;
        }
        const el = document.getElementById('phase-plot');
        if (!el) {
            console.error('plotPhaseShifts: div #phase-plot not found');
            return;
        }
        const traces = {};
        data.forEach(point => {
            const L = point.L;
            const phaseShift = point.phase_shift != null ? point.phase_shift : point.phaseShift;
            const energy = point.energy;
            if (L === undefined || (phaseShift === undefined && phaseShift !== 0)) return;
            if (!traces[L]) {
                traces[L] = {
                    x: [],
                    y: [],
                    name: `L = ${L}`,
                    type: 'scatter',
                    mode: 'lines+markers',
                    line: { width: 3, ...PLOTLY_SPLINE_LINE },
                    marker: { size: 6 }
                };
            }
            traces[L].x.push(Number(energy));
            traces[L].y.push((Number(phaseShift) || 0) * 180 / Math.PI);
        });

        const plotData = Object.values(traces);
        plotData.forEach(sortTraceByX);
        if (plotData.length === 0) {
            console.warn('plotPhaseShifts: no traces after processing');
            return;
        }
        const layout = {
            title: 'Nuclear Phase Shifts vs Energy',
            xaxis: { title: 'Energy (MeV)', gridcolor: '#e0e0e0' },
            yaxis: { title: 'Phase Shift (degrees)', gridcolor: '#e0e0e0' },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Arial, sans-serif' },
            legend: { x: 0.02, y: 0.98 },
            margin: { t: 50, b: 50, l: 60, r: 30 }
        };

        try {
            if (typeof Plotly === 'undefined') throw new Error('Plotly not loaded');
            Plotly.newPlot('phase-plot', plotData, layout, {responsive: true});
        } catch (err) {
            console.error('Plotly.newPlot failed:', err);
        }
    }

    plotRMatrices() {
        const data = this.currentData.r_matrices || this.currentData['r-matrices'];
        if (!data || !Array.isArray(data) || data.length === 0) return;
        const traces = {};

        data.forEach(point => {
            const L = point.L;
            if (!traces[L]) {
                traces[L] = {
                    nuclear: { x: [], y: [], name: `L = ${L} (Nuclear)`, type: 'scatter', mode: 'lines+markers', line: { ...PLOTLY_SPLINE_LINE } },
                    coulomb_nuclear: { x: [], y: [], name: `L = ${L} (Coul+Nuc)`, type: 'scatter', mode: 'lines+markers', line: { ...PLOTLY_SPLINE_LINE, dash: 'dash' } }
                };
            }
            traces[L].nuclear.x.push(point.energy);
            traces[L].nuclear.y.push(point.r_nuclear);
            traces[L].coulomb_nuclear.x.push(point.energy);
            traces[L].coulomb_nuclear.y.push(point.r_coulomb_nuclear);
        });

        const plotData = [];
        Object.values(traces).forEach(trace => {
            sortTraceByX(trace.nuclear);
            sortTraceByX(trace.coulomb_nuclear);
            plotData.push(trace.nuclear, trace.coulomb_nuclear);
        });

        const layout = {
            title: 'R-Matrix Values Comparison',
            xaxis: { title: 'Energy (MeV)', gridcolor: '#e0e0e0' },
            yaxis: { title: 'R-Matrix', gridcolor: '#e0e0e0' },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Arial, sans-serif' },
            legend: { x: 0.02, y: 0.98 },
            margin: { t: 50, b: 50, l: 60, r: 30 }
        };

        Plotly.newPlot('rmatrix-plot', plotData, layout, {responsive: true});
    }

    plotPotentials() {
        const data = this.currentData.potentials;
        if (!data || !Array.isArray(data) || data.length === 0) return;

        const woodsSaxon = {
            x: data.map(p => p.radius),
            y: data.map(p => p.woods_saxon),
            name: 'Woods-Saxon',
            type: 'scatter',
            mode: 'lines',
            line: { color: 'blue', width: 3, ...PLOTLY_SPLINE_LINE }
        };

        const coulomb = {
            x: data.map(p => p.radius),
            y: data.map(p => p.coulomb),
            name: 'Coulomb',
            type: 'scatter',
            mode: 'lines',
            line: { color: 'red', width: 3, ...PLOTLY_SPLINE_LINE }
        };

        const combined = {
            x: data.map(p => p.radius),
            y: data.map(p => p.combined),
            name: 'Combined',
            type: 'scatter',
            mode: 'lines',
            line: { color: 'green', width: 3, ...PLOTLY_SPLINE_LINE }
        };

        const layout = {
            title: 'Nuclear Potentials vs Radius',
            xaxis: { title: 'Radius (fm)', gridcolor: '#e0e0e0' },
            yaxis: { title: 'Potential (MeV)', gridcolor: '#e0e0e0' },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Arial, sans-serif' },
            legend: { x: 0.02, y: 0.98 },
            margin: { t: 50, b: 50, l: 60, r: 30 }
        };

        Plotly.newPlot('potential-plot', [woodsSaxon, coulomb, combined], layout, {responsive: true});
    }

    plotCrossSections() {
        const data = this.currentData.cross_sections || this.currentData['cross-sections'];
        if (!data || !Array.isArray(data) || data.length === 0) return;

        const trace = {
            x: data.map(p => p.energy),
            y: data.map(p => p.total_cross_section),
            name: 'Total Cross-Section',
            type: 'scatter',
            mode: 'lines+markers',
            line: { color: 'purple', width: 3, ...PLOTLY_SPLINE_LINE },
            marker: { size: 6 }
        };
        sortTraceByX(trace);

        const layout = {
            title: 'Total Cross-Sections vs Energy',
            xaxis: { title: 'Energy (MeV)', gridcolor: '#e0e0e0' },
            yaxis: { title: 'Cross-Section (arbitrary units)', type: 'log', gridcolor: '#e0e0e0' },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Arial, sans-serif' },
            margin: { t: 50, b: 50, l: 60, r: 30 }
        };

        Plotly.newPlot('cross-section-plot', [trace], layout, {responsive: true});
    }

    plotDashboard() {
        // Create a comprehensive dashboard with multiple subplots
        const phaseData = this.currentData.phase_shifts || this.currentData['phase-shifts'];
        const potentialData = this.currentData.potentials;
        const crossSectionData = this.currentData.cross_sections || this.currentData['cross-sections'];
        if (!phaseData?.length || !potentialData?.length || !crossSectionData?.length) return;

        // Phase shifts (grouped by L)
        const phaseTraces = {};
        phaseData.forEach(point => {
            const L = point.L;
            if (!phaseTraces[L]) {
                phaseTraces[L] = {
                    x: [],
                    y: [],
                    name: `L = ${L}`,
                    type: 'scatter',
                    mode: 'lines+markers',
                    line: { ...PLOTLY_SPLINE_LINE },
                    showlegend: true
                };
            }
            phaseTraces[L].x.push(point.energy);
            phaseTraces[L].y.push(point.phase_shift * 180 / Math.PI);
        });

        Object.values(phaseTraces).forEach(sortTraceByX);

        const dashWoods = {
            x: potentialData.map(p => p.radius),
            y: potentialData.map(p => p.woods_saxon),
            name: 'Woods-Saxon',
            type: 'scatter',
            mode: 'lines',
            line: { ...PLOTLY_SPLINE_LINE },
            xaxis: 'x2',
            yaxis: 'y2',
            showlegend: false
        };
        const dashCoul = {
            x: potentialData.map(p => p.radius),
            y: potentialData.map(p => p.coulomb),
            name: 'Coulomb',
            type: 'scatter',
            mode: 'lines',
            line: { ...PLOTLY_SPLINE_LINE },
            xaxis: 'x2',
            yaxis: 'y2',
            showlegend: false
        };
        const dashXs = {
            x: crossSectionData.map(p => p.energy),
            y: crossSectionData.map(p => p.total_cross_section),
            name: 'Cross-Section',
            type: 'scatter',
            mode: 'lines',
            line: { ...PLOTLY_SPLINE_LINE },
            xaxis: 'x3',
            yaxis: 'y3',
            showlegend: false
        };
        sortTraceByX(dashWoods);
        sortTraceByX(dashCoul);
        sortTraceByX(dashXs);

        const traces = [
            // Phase shifts
            ...Object.values(phaseTraces),
            dashWoods,
            dashCoul,
            dashXs
        ];

        const layout = {
            title: 'DWBA Comprehensive Dashboard',
            grid: {
                rows: 2,
                columns: 2,
                subplots: [
                    ['xy', 'x2y2'],
                    ['xy', 'x3y3']
                ]
            },
            xaxis: { title: 'Energy (MeV)', domain: [0, 0.45] },
            yaxis: { title: 'Phase Shift (degrees)', domain: [0.55, 1] },
            xaxis2: { title: 'Radius (fm)', domain: [0.55, 1] },
            yaxis2: { title: 'Potential (MeV)', domain: [0.55, 1] },
            xaxis3: { title: 'Energy (MeV)', domain: [0, 0.45] },
            yaxis3: { title: 'Cross-Section', type: 'log', domain: [0, 0.45] },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Arial, sans-serif' },
            margin: { t: 50, b: 30, l: 50, r: 30 }
        };

        Plotly.newPlot('dashboard-plot', traces, layout, {responsive: true});
    }

    plotElastic() {
        const data = this.currentData.elastic;
        if (!data || data.length === 0) return;

        const ratioMode = !!(document.getElementById('elastic_ratio_rutherford') || {}).checked;
        const projK = this.elasticProjectileKinematics();
        const targAZ = this.elasticTargetAZ();
        const mTarg = targAZ.A * 931.5;
        const z1 = projK.Z;
        const z2 = targAZ.Z;
        const canRatio = ratioMode && z1 * z2 > 0;
        const ratioUnavailable = ratioMode && z1 * z2 <= 0;

        // Get current lab energies from Elastic tab to filter data
        const energyEl = document.getElementById('elastic-lab-energies');
        const requestedEnergies = energyEl && energyEl.value
            ? energyEl.value.split(',').map(e => parseFloat(e.trim())).filter(e => !isNaN(e))
            : null;

        const traces = {};
        
        // Group by energy, optionally filtering to requested energies (float-tolerant; JSON may differ slightly)
        const energyRequestedP = (e) => {
            if (!requestedEnergies) return true;
            const x = Number(e);
            return requestedEnergies.some((re) => Math.abs(re - x) <= 1e-5 * Math.max(1, Math.abs(re)) || Math.abs(re - x) < 1e-4);
        };
        data.forEach(point => {
            const E = point.energy;
            if (!energyRequestedP(E)) return;
            if (!traces[E]) {
                traces[E] = {
                    x: [],
                    y: [],
                    name: `E = ${E} MeV`,
                    type: 'scatter',
                    mode: 'lines+markers',
                    line: { width: 2, ...PLOTLY_SPLINE_LINE },
                    marker: { size: 4 }
                };
            }
            const dsigma = point.differential_cross_section;
            let yVal = dsigma;
            if (canRatio) {
                const eCm = this.labToCmKineticMeV(projK.mass, mTarg, E);
                const ruth = this.rutherfordDsDomegaMbSr(z1, z2, eCm, point.angle);
                yVal = ruth && ruth > 0 ? dsigma / ruth : NaN;
            }
            traces[E].x.push(point.angle);
            traces[E].y.push(yVal);
        });

        const plotData = Object.values(traces);
        plotData.forEach(sortTraceByX);
        const layout = {
            title: canRatio
                ? 'Elastic dσ/dΩ / σ_Rutherford'
                : ratioUnavailable
                    ? 'Elastic dσ/dΩ (mb/sr)'
                    : 'Elastic differential cross section',
            xaxis: { title: 'Scattering angle θ_cm (degrees)', gridcolor: '#e0e0e0' },
            yaxis: canRatio
                ? { title: 'dσ/dΩ / σ_Rutherford', type: 'log', gridcolor: '#e0e0e0' }
                : { title: 'dσ/dΩ (mb/sr)', type: 'log', gridcolor: '#e0e0e0' },
            annotations: [],
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Arial, sans-serif' },
            legend: { x: 0.02, y: 0.98 },
            margin: { t: 50, b: 50, l: 60, r: 30 }
        };

        Plotly.newPlot('elastic-plot', plotData, layout, {responsive: true});
    }

    plotInelastic() {
        const data = this.currentData.inelastic;
        if (!data || data.length === 0) return;

        const hasAngle = data[0] && data[0].angle !== undefined;
        const logFloor = 1e-40;

        if (hasAngle) {
            const impl = data[0].implementation;
            const energies = [...new Set(data.map(p => p.energy))].sort((a, b) => a - b);
            const plotData = energies.map(E => {
                const pts = data.filter(p => p.energy === E);
                const tr = {
                    x: pts.map(p => p.angle),
                    y: pts.map(p => Math.max(p.differential_cross_section, logFloor)),
                    name: energies.length > 1 ? `E_p,lab = ${E} MeV` : `E_p,lab = ${E} MeV`,
                    type: 'scatter',
                    mode: 'lines+markers',
                    line: { width: 3, ...PLOTLY_SPLINE_LINE },
                    marker: { size: 6 }
                };
                sortTraceByX(tr);
                return tr;
            });
            const paperTitle = impl === 'li11_paper1708' ? '¹¹Li (p,p′) paper1708 — ' : '';
            const layout = {
                title: paperTitle + 'Inelastic dσ/dΩ (CM, mb/sr)',
                xaxis: { title: 'Scattering angle θ_cm (degrees)', gridcolor: '#e0e0e0' },
                yaxis: { title: 'dσ/dΩ (mb/sr)', type: 'log', gridcolor: '#e0e0e0' },
                plot_bgcolor: 'rgba(0,0,0,0)',
                paper_bgcolor: 'rgba(0,0,0,0)',
                font: { family: 'Arial, sans-serif' },
                legend: { x: 0.02, y: 0.98 },
                margin: { t: 50, b: 50, l: 60, r: 30 }
            };
            Plotly.newPlot('inelastic-plot', plotData, layout, { responsive: true });
            return;
        }

        const traces = {};
        const groupByLambda = data[0] && data[0].lambda !== undefined;

        data.forEach(point => {
            const key = groupByLambda ? point.lambda : point.L;
            const label = groupByLambda ? `λ = ${key}` : `L = ${key}`;
            if (!traces[key]) {
                traces[key] = {
                    x: [],
                    y: [],
                    name: label,
                    type: 'scatter',
                    mode: 'lines+markers',
                    line: { width: 3, ...PLOTLY_SPLINE_LINE },
                    marker: { size: 6 }
                };
            }
            traces[key].x.push(point.energy);
            const y = Number(point.differential_cross_section);
            traces[key].y.push(Number.isFinite(y) ? Math.max(y, logFloor) : logFloor);
        });

        const plotData = Object.values(traces);
        plotData.forEach(sortTraceByX);
        const projEl = document.getElementById('inelastic_projectile');
        const tgtEl = document.getElementById('inelastic_target');
        const projLabels = { p: 'p', n: 'n', d: 'd', a: 'α' };
        const tgtLabels = { '12C': '¹²C', '16O': '¹⁶O', generic: 'target' };
        const proj = projLabels[projEl?.value] || projEl?.value || 'p';
        const tgt = tgtLabels[tgtEl?.value] || tgtEl?.value || '';
        const reactionLabel = tgt ? `${proj} + ${tgt} ` : '';
        const layout = {
            title: reactionLabel + (groupByLambda ? 'inelastic dσ/dΩ by multipole λ (summed over L)' : 'Inelastic Scattering Differential Cross-Section'),
            xaxis: { title: 'Energy (MeV)', gridcolor: '#e0e0e0' },
            yaxis: { title: 'dσ/dΩ (mb/sr)', type: 'log', gridcolor: '#e0e0e0' },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Arial, sans-serif' },
            legend: { x: 0.02, y: 0.98 },
            margin: { t: 50, b: 50, l: 60, r: 30 }
        };

        Plotly.newPlot('inelastic-plot', plotData, layout, {responsive: true});
    }

    plotTransfer() {
        const data = this.currentData.transfer;
        if (!data || data.length === 0) return;

        // Data can be DCS vs. angle (has angle) or legacy DCS vs. E per L
        const hasAngle = data[0] && data[0].angle !== undefined;

        let plotData, xTitle;
        if (hasAngle) {
            xTitle = 'Scattering angle (CM, degrees)';
            const energies = [...new Set(data.map(p => p.energy))].sort((a, b) => a - b);
            const label = this.currentData.transferTargetLabel || 'Transfer';
            // Log scale can't show 0; floor tiny values so L=1 node at 90° doesn't drop off chart
            const logFloor = 1e-40;
            plotData = energies.map(E => {
                const pts = data.filter(p => p.energy === E);
                const tr = {
                    x: pts.map(p => p.angle),
                    y: pts.map(p => Math.max(p.differential_cross_section, logFloor)),
                    name: energies.length > 1 ? `E = ${E} MeV` : `${label}, E = ${E} MeV`,
                    type: 'scatter',
                    mode: 'lines+markers',
                    line: { width: 3, ...PLOTLY_SPLINE_LINE },
                    marker: { size: 6 }
                };
                sortTraceByX(tr);
                return tr;
            });
        } else {
            const traces = {};
            data.forEach(point => {
                const L = point.L;
                if (!traces[L]) {
                    traces[L] = {
                        x: [],
                        y: [],
                        name: `L = ${L}`,
                        type: 'scatter',
                        mode: 'lines+markers',
                        line: { width: 3, ...PLOTLY_SPLINE_LINE },
                        marker: { size: 6 }
                    };
                }
                traces[L].x.push(point.energy);
                traces[L].y.push(point.differential_cross_section);
            });
            plotData = Object.values(traces);
            plotData.forEach(sortTraceByX);
            xTitle = 'Energy (MeV)';
        }

        const layout = {
            title: 'Transfer Reaction Differential Cross-Section',
            xaxis: { title: xTitle, gridcolor: '#e0e0e0' },
            yaxis: { title: 'dσ/dΩ (mb/sr)', type: 'log', gridcolor: '#e0e0e0' },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Arial, sans-serif' },
            legend: { x: 0.02, y: 0.98 },
            margin: { t: 50, b: 50, l: 60, r: 30 }
        };

        Plotly.newPlot('transfer-plot', plotData, layout, {responsive: true});
    }

    updateDashboardStats(calculationTime) {
        if (!this.currentData) return;
        const phaseShifts = this.currentData.phase_shifts || this.currentData['phase-shifts'];
        if (!phaseShifts?.length) return;

        const totalPoints = phaseShifts.length;
        const energies = phaseShifts.map(p => p.energy);
        const LValues = [...new Set(phaseShifts.map(p => p.L))];

        document.getElementById('total-points').textContent = totalPoints;
        document.getElementById('energy-range-display').textContent = 
            `${Math.min(...energies)}-${Math.max(...energies)}`;
        document.getElementById('L-count').textContent = LValues.length;
        document.getElementById('calculation-time').textContent = `${calculationTime}ms`;
    }

    resetParameters() {
        this.loadDefaultParameters();
        this.showStatus('Parameters reset to defaults', 'info');
    }
}

let dashboard;
document.addEventListener('DOMContentLoaded', () => {
    dashboard = new DWBADashboard();
});
