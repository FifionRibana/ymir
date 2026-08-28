//! Thread-local cancel token for cooperative mid-step interruption.
//!
//! The v2 bridge (in `crates/ymir-viz/src/bridge/v2`) owns an
//! `Arc<AtomicBool>` that the UI flips to `true` when the user presses
//! the Stop button. Before Step 12 Phase 7b, that flag was only
//! observed at the *outer* step boundary of
//! [`crate::tectonics_v2::diagnostics::harness::run_baseline_with_progress`]
//! — pressing Stop on a 64² mantle-on run left the user waiting one
//! step (≈ 5–25 s) for the current CG / Newton iter to finish.
//!
//! This module exposes the same flag as a thread-local for the inner
//! loops of the solver and the HD erosion to consult. The deployment
//! model is "one solver thread per bridge", so a thread-local is
//! sufficient and avoids adding an `Option<Arc<AtomicBool>>` argument
//! through every solver call signature.
//!
//! Lifetime:
//!
//! - At the start of each `V2Command` handler in
//!   `crates/ymir-viz/src/bridge/v2/thread.rs`, the bridge calls
//!   [`set`] with `Some(cancel.clone())`. The cancel flag is now
//!   observable from any code reached on this thread.
//! - During execution, the CG inner loop, the Newton outer loop, the
//!   harness step loop, and the HD erosion callback wrapper all call
//!   [`is_cancelled`]. None of them mutate the token.
//! - At command end, the handler calls [`clear`] (or [`set`] with
//!   `None`) to drop the token so a future independent run does not
//!   inherit a stale cancel signal.
//!
//! When no token is bound, [`is_cancelled`] returns `false` (the
//! "running outside the bridge" semantics — e.g., unit tests calling
//! `run_baseline` directly). Bit-determinism: a `false`-checked token
//! never changes the trajectory, so a run with no token bound and a
//! run with a token bound that stays `false` produce byte-identical
//! `BaselineResult`s.

use std::cell::RefCell;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

thread_local! {
    static CURRENT_TOKEN: RefCell<Option<Arc<AtomicBool>>> = const { RefCell::new(None) };
}

/// Bind a cancel token to the current thread. Pass `None` to clear.
/// The caller is responsible for clearing at the end of a logical
/// run boundary; subsequent runs on the same thread that forget to
/// clear would observe the prior run's flag.
pub fn set(token: Option<Arc<AtomicBool>>) {
    CURRENT_TOKEN.with(|cell| *cell.borrow_mut() = token);
}

/// Convenience: clear the bound token. Equivalent to `set(None)`.
pub fn clear() {
    set(None);
}

/// Returns `true` iff a cancel token is currently bound on this
/// thread *and* the underlying [`AtomicBool`] is set.
///
/// Reads use `Ordering::Relaxed`: cancel is advisory, not a
/// memory-ordering primitive, and inner loops call this frequently
/// (once per ~16 CG iters, once per Newton outer iter) so the
/// cheapest load is preferred. The bridge sets the flag via
/// `Ordering::Relaxed` for the same reason; the consequent worst-
/// case delay between a Stop click and the next observed `true` is
/// a few microseconds, dominated by the inner-loop iteration
/// cadence rather than the atomic ordering.
pub fn is_cancelled() -> bool {
    CURRENT_TOKEN
        .with(|cell| cell.borrow().as_ref().map(|a| a.load(Ordering::Relaxed)).unwrap_or(false))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// No token bound → `is_cancelled` returns `false`.
    #[test]
    fn no_token_returns_false() {
        clear();
        assert!(!is_cancelled());
    }

    /// Token bound but flag false → `is_cancelled` returns `false`.
    #[test]
    fn bound_token_with_false_flag_returns_false() {
        let token = Arc::new(AtomicBool::new(false));
        set(Some(token.clone()));
        assert!(!is_cancelled());
        clear();
    }

    /// Token bound and flag set externally → `is_cancelled` returns
    /// `true`. Mirrors the bridge → solver-thread interaction:
    /// `set` is called by the solver-thread side at command entry,
    /// `store(true)` is called by the UI-thread side on Stop click.
    #[test]
    fn external_flip_observable() {
        let token = Arc::new(AtomicBool::new(false));
        set(Some(token.clone()));
        token.store(true, Ordering::Relaxed);
        assert!(is_cancelled());
        clear();
    }

    /// `clear` drops the bound token. A subsequent `is_cancelled`
    /// returns `false` even if the externally-held `Arc` still has
    /// `true` set. Mirrors the bridge's "end of command" cleanup.
    #[test]
    fn clear_drops_observation() {
        let token = Arc::new(AtomicBool::new(true));
        set(Some(token.clone()));
        assert!(is_cancelled());
        clear();
        assert!(!is_cancelled());
        // The externally-held Arc still sees true; it just isn't
        // observed here anymore.
        assert!(token.load(Ordering::Relaxed));
    }
}
