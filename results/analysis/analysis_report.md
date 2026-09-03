# Verified presented-paper analysis

- Bennu/Ryugu full rerun triggered: False
- Literature-validation cases: 18

## Standalone literature validation

| Body | Method | Metric | Reference | N | MAE (m) | MAPE (%) |
|---|---|---|---|---:|---:|---:|
| Bennu | APBT | depth | Daly | 33 | 1.343 | 28.938 |
| Bennu | APBT | diameter | Daly | 33 | 6.504 | 10.757 |
| Bennu | APBT | diameter | Bierhaus | 44 | 6.251 | 13.266 |
| Bennu | APBT | diameter | Deshapriya | 45 | 6.749 | 13.387 |
| Ryugu | APBT | depth | Noguchi | 77 | 0.916 | 17.976 |
| Ryugu | APBT | diameter | Noguchi | 77 | 5.551 | 11.464 |
| Itokawa | APBT | depth | Naru-Hirata | 21 | 2.663 | 66.688 |
| Itokawa | APBT | diameter | Naru-Hirata | 21 | 5.104 | 9.832 |
| Didymos | APBT | diameter | Barnouin | 16 | 11.597 | 12.215 |
| Bennu | CIRCLE | depth | Daly | 32 | 1.948 | 38.803 |
| Bennu | CIRCLE | diameter | Daly | 32 | 8.003 | 14.051 |
| Bennu | CIRCLE | diameter | Bierhaus | 41 | 7.001 | 18.278 |
| Bennu | CIRCLE | diameter | Deshapriya | 42 | 8.115 | 19.222 |
| Ryugu | CIRCLE | depth | Noguchi | 77 | 1.123 | 21.577 |
| Ryugu | CIRCLE | diameter | Noguchi | 77 | 6.747 | 14.228 |
| Itokawa | CIRCLE | depth | Naru-Hirata | 20 | 4.143 | 82.063 |
| Itokawa | CIRCLE | diameter | Naru-Hirata | 20 | 7.527 | 13.257 |
| Didymos | CIRCLE | diameter | Barnouin | 14 | 17.854 | 16.102 |

## Source provenance

- Bennu APBT: `results/Bennu/apbt.json` (declared_final_release).
- Bennu CIRCLE: `results/Bennu/circle.json` (declared_final_release).
- Ryugu APBT: `results/Ryugu/apbt.json` (declared_final_release).
- Ryugu CIRCLE: `results/Ryugu/circle.json` (declared_final_release).
- Itokawa APBT: `results/Itokawa/apbt.json` (declared_final_release).
- Itokawa CIRCLE: `results/Itokawa/circle.json` (declared_final_release).
- Didymos APBT: `results/Didymos/apbt.json` (declared_final_release).
- Didymos CIRCLE: `results/Didymos/circle.json` (declared_final_release).
