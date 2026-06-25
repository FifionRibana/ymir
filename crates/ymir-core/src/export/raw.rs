//! Low-level raw binary codec — `f32` / `u8` / `u32`, little-endian, row-major,
//! **no header**.
//!
//! This is the SINGLE binary codec, SHARED by two consumers with distinct
//! policies (see `docs/design/c1_generation_cache.md`):
//!   - the **§9 export** (`PipelineExport`) — a stable, human-named Living Landz
//!     deliverable whose public contract is `metadata.json` (§9.3);
//!   - the **generation cache** (`crate::cache`) — internal, disposable,
//!     content-addressed scratch for skipping unchanged upstream steps.
//!
//! Both write the exact same bytes (`val.to_le_bytes()` concatenated, no
//! padding, no magic). Keeping ONE codec here prevents a second serializer from
//! drifting away from the §9 byte layout. `GridF32::save_raw` / `load_raw`
//! delegate here too, so every `.raw` in the project is the same format.

use std::path::Path;

/// Write `f32` slice as little-endian bytes (no header).
pub fn save_f32(path: &Path, data: &[f32]) -> Result<(), String> {
    let bytes: Vec<u8> = data.iter().flat_map(|f| f.to_le_bytes()).collect();
    std::fs::write(path, bytes).map_err(|e| format!("Write error: {e}"))
}

/// Read `expected_len` little-endian `f32` values (errors on size mismatch).
pub fn load_f32(path: &Path, expected_len: usize) -> Result<Vec<f32>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("Read error: {e}"))?;
    if bytes.len() != expected_len * 4 {
        return Err(format!(
            "Size mismatch: expected {} bytes, got {}",
            expected_len * 4,
            bytes.len()
        ));
    }
    Ok(bytes.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect())
}

/// Write `u8` slice (identity — bytes are already little-endian).
pub fn save_u8(path: &Path, data: &[u8]) -> Result<(), String> {
    std::fs::write(path, data).map_err(|e| format!("Write error: {e}"))
}

/// Read all bytes as `u8`.
pub fn load_u8(path: &Path) -> Result<Vec<u8>, String> {
    std::fs::read(path).map_err(|e| format!("Read error: {e}"))
}

/// Write `u32` slice as little-endian bytes (no header).
pub fn save_u32(path: &Path, data: &[u32]) -> Result<(), String> {
    let bytes: Vec<u8> = data.iter().flat_map(|v| v.to_le_bytes()).collect();
    std::fs::write(path, bytes).map_err(|e| format!("Write error: {e}"))
}

/// Read `expected_len` little-endian `u32` values (errors on size mismatch).
pub fn load_u32(path: &Path, expected_len: usize) -> Result<Vec<u32>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("Read error: {e}"))?;
    if bytes.len() != expected_len * 4 {
        return Err(format!(
            "Size mismatch: expected {} bytes, got {}",
            expected_len * 4,
            bytes.len()
        ));
    }
    Ok(bytes.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f32_round_trip_exact() {
        let dir = std::env::temp_dir().join("ymir_raw_codec_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("rt.raw");
        let data = vec![0.0f32, -1.5, 3.1415927, f32::MIN_POSITIVE, 12345.678];
        save_f32(&path, &data).unwrap();
        let back = load_f32(&path, data.len()).unwrap();
        assert_eq!(data, back, "f32 raw round-trip must be byte-exact");
        // wrong length → error, never a silent truncation.
        assert!(load_f32(&path, data.len() + 1).is_err());
    }

    #[test]
    fn u32_round_trip_exact() {
        let dir = std::env::temp_dir().join("ymir_raw_codec_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("rt_u32.raw");
        let data = vec![0u32, 1, 4_294_967_295, 42];
        save_u32(&path, &data).unwrap();
        assert_eq!(data, load_u32(&path, data.len()).unwrap());
    }
}
