# Verified presented-paper analysis

- Bennu/Ryugu full rerun triggered: False
- Literature-validation cases: 20

## Standalone literature validation

| Body | Method | Metric | Reference | N | MAE (m) | MAPE (%) |
|---|---|---|---|---:|---:|---:|
| Bennu | APBT | depth | Daly | 33 | 1.340 | 28.886 |
| Bennu | APBT | diameter | Daly | 33 | 6.590 | 10.842 |
| Bennu | APBT | diameter | Bierhaus | 44 | 6.186 | 13.210 |
| Bennu | APBT | diameter | Deshapriya | 45 | 6.686 | 13.332 |
| Ryugu | APBT | depth | Noguchi | 77 | 0.916 | 17.978 |
| Ryugu | APBT | diameter | Noguchi | 77 | 5.569 | 11.482 |
| Ryugu | APBT | diameter | Hirata | 77 | 7.568 | 22.306 |
| Itokawa | APBT | depth | Naru-Hirata | 21 | 2.690 | 66.606 |
| Itokawa | APBT | diameter | Naru-Hirata | 21 | 4.976 | 9.657 |
| Didymos | APBT | diameter | Barnouin | 16 | 11.682 | 12.247 |
| Bennu | CIRCLE | depth | Daly | 32 | 1.966 | 39.058 |
| Bennu | CIRCLE | diameter | Daly | 32 | 8.128 | 14.175 |
| Bennu | CIRCLE | diameter | Bierhaus | 41 | 6.903 | 18.193 |
| Bennu | CIRCLE | diameter | Deshapriya | 42 | 8.019 | 19.139 |
| Ryugu | CIRCLE | depth | Noguchi | 77 | 1.123 | 21.574 |
| Ryugu | CIRCLE | diameter | Noguchi | 77 | 6.746 | 14.226 |
| Ryugu | CIRCLE | diameter | Hirata | 77 | 9.294 | 30.813 |
| Itokawa | CIRCLE | depth | Naru-Hirata | 19 | 3.719 | 83.562 |
| Itokawa | CIRCLE | diameter | Naru-Hirata | 19 | 7.301 | 13.458 |
| Didymos | CIRCLE | diameter | Barnouin | 14 | 17.531 | 15.982 |

## Source provenance

- Bennu APBT: `results/Bennu/apbt.json` (declared_final_release).
- Bennu CIRCLE: `results/Bennu/circle.json` (declared_final_release).
- Ryugu APBT: `results/Ryugu/apbt.json` (declared_final_release).
- Ryugu CIRCLE: `results/Ryugu/circle.json` (declared_final_release).
- Itokawa APBT: `results/Itokawa/apbt.json` (declared_final_release).
- Itokawa CIRCLE: `results/Itokawa/circle.json` (declared_final_release).
- Didymos APBT: `results/Didymos/apbt.json` (declared_final_release).
- Didymos CIRCLE: `results/Didymos/circle.json` (declared_final_release).
