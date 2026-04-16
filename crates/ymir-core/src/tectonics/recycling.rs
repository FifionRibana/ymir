//! Conservative mass recycling: subducted material is redistributed
//! to volcanic arcs (immediately) and spreading ridges (after delay).
//!
//! All creation at convergent/divergent boundaries is funded by
//! subduction destruction — total mass is conserved (minus optional loss).

/// Configuration for the mass recycling system.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RecyclingConfig {
    /// Whether conservative recycling is enabled.
    /// When enabled, spreading and volcanic arc rates are ignored —
    /// all creation is funded by the recycling budget.
    /// Default: true.
    pub enabled: bool,
    /// Fraction of subducted mass that goes immediately to volcanic arc
    /// on the overriding plate. Models fast magma ascent from the slab.
    /// Default: 0.15.
    pub arc_fraction: f64,
    /// Fraction of subducted mass permanently lost to the deep mantle.
    /// Default: 0.0 (full recycling).
    pub loss_fraction: f64,
    /// Delay in timesteps before subducted mass becomes available for
    /// spreading. Models mantle transit time. Default: 20.
    pub mantle_delay: usize,
}

impl Default for RecyclingConfig {
    fn default() -> Self {
        Self { enabled: true, arc_fraction: 0.15, loss_fraction: 0.0, mantle_delay: 20 }
    }
}

/// Tracks subducted mass and distributes it over time.
pub struct RecyclingBuffer {
    buffer: Vec<f64>,
    head: usize,
    capacity: usize,
}

impl RecyclingBuffer {
    pub fn new(delay: usize) -> Self {
        let capacity = delay.max(1);
        Self { buffer: vec![0.0; capacity], head: 0, capacity }
    }

    /// Deposit mass into the buffer at the current head position.
    /// It will become available after `delay` steps (when the head
    /// wraps back to this slot).
    pub fn deposit(&mut self, mass: f64) {
        let write_pos = self.head % self.capacity;
        self.buffer[write_pos] += mass;
    }

    /// Advance one step and return the mass that becomes available
    /// for spreading. The returned slot is cleared.
    pub fn advance(&mut self) -> f64 {
        let read_pos = (self.head + 1) % self.capacity;
        let available = self.buffer[read_pos];
        self.buffer[read_pos] = 0.0;
        self.head = read_pos;
        available
    }
}
