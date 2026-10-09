# F163 — predictions, written before any measurement (2026-10-09)

**Not blind**:
- F161's P0 readings at 256² (λ 17.64 km, facets 13.0 %, octaves 32.6 / 121.4 m against Corsica's 84.9 / 152.4 m);
- F162's 2 048² readings and images;
- the code facts of `f163_declared.md` D0 (the MFD area is already on; the head threshold and the diffusion are
  sub-cell).

Nothing about the three new instruments, the variants or the controls has been run.

## The reviewer's (from the brief, verbatim in substance)

1. On P0 the directions are ≥ 2× more aligned than Corsica's, and the trunk sinuosity is ≥ 20 % lower.
2. P-mfd (or the random receiver) brings the directions within ×1.5 of Corsica.
3. P-dissection brings λ and the 3–6 km octave closer and lowers the physics level's facets.
4. The best combination through R4+π misses no octave beyond ×1.5 at 2 048².
5. Meta: at least one of these is false.

**How they are scored, declared**:
- (1) on P0 at 256² and on P0's chain at 2 048²; the sinuosity on S_W itself, as written, and on S_W − 1 too;
- (2) and (3) at 256²;
- (4) on the 0.4–12.5 km octaves at 2 048².

## Mine

**Controls**:
- **C1**: the D8 plane at 10° reads A_dir 45–70 %; Corsica 25–35 % at every cell; the fractal 22–30 %. The
  instrument passes, though not by much.
- **C2**: the straight valley reads S − 1 ≤ 0.01 at 0°, and 0.03–0.08 at 22.5° (the D8 floor).
  - Corsica's S_5 − 1 at 195 m is 0.15–0.35. The instrument passes at 195 m.
  - At 1 563 m (S_10 only) it **fails** its control: the D8 floor is too close to Corsica's coarse sinuosity, so it is
    blind there.
- **C3**: the terraced field reads F ≥ 60 %; Corsica 5–15 % at 1 563 m and 1–5 % at 195 m. The instrument passes.

**Measures**:
- **M1, P0 at 256²**: A_dir 40–55 % (1.3–1.8 × Corsica); S_10 − 1 of 0.02–0.06; F 20–35 % (outside ×1.5 of Corsica).
- **M2, N1 and R4+π at 2 048²**:
  - A_dir 40–55 % against Corsica's 25–35 %;
  - S_5 − 1 of 0.05–0.10 against 0.15–0.30;
  - F 5–15 % against 1–5 %: outside ×1.5 on all three.

**The variants at 256²**:
- **D1, P-mfd (the random receiver)**:
  - A_dir within ×1.5 of Corsica;
  - λ within ±15 % of P0's;
  - the facets down by ≥ 20 %;
  - each step 20–60 % dearer.
- **D2, P-bruit**:
  - λ down by 10–25 %;
  - the 3.1–6.25 km octave up by 10–40 %;
  - share (a) 5–15 %.
- **D3, P-dissection**: no p brings λ within ×1.25 of 9.08 km; D8 (no MFD) gives facets ≥ P0's.

**D4, at 2 048²**:
- the best variant is a combination that includes the random receiver;
- through R4+π it is **rejected at 2 048²**:
  - on the 1.5–6 km octaves (F162's diagnosis unchanged);
  - and on F;
- A_dir at 2 048² closes less than half of the gap it closes at 256², because the fine levels carve their own D8
  streams.

**C, the cost**:
- each physics variant takes 5–40 s;
- the two R4+π chains take about 4 min each;
- the whole bench takes under 30 min.

**Meta**: at least two of mine are refuted.
