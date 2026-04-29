This directory is reserved for scan tables, appended SuperIso outputs, and run logs.

The reusable helper programs that operate on these data files are kept outside `np_data/` and now live in `../third_party_patches/superiso_v5.0/`:

- `../third_party_patches/superiso_v5.0/check_superiso.py`
- `../third_party_patches/superiso_v5.0/run_superiso_batch.sh`

Typical calls from the repository root are:

```bash
python3 -u third_party_patches/superiso_v5.0/check_superiso.py np_data/scan_after_hs.dat \
	--mode chi2 --append-table --jobs 8 --checkpoint-every 50 \
	--output-file np_data/scan_after_hs_superiso.dat

bash third_party_patches/superiso_v5.0/run_superiso_batch.sh 8 50
```
