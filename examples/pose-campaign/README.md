# Uncomputed pose campaign

This is a ready-to-run set of nine explicit 53-atom starting candidates: axial separations 3.4, 3.6 and 3.8 Å, crossed with lateral offsets −0.2, 0 and +0.2 Å. All use PBE0-D3(BJ)/def2-SVP, a neutral doublet, and the same six distal anchors. No calculation has been run in this example campaign.

Copy this directory to a new working location before executing it:

```bash
mkdir -p campaigns
cp -R examples/pose-campaign campaigns/my-poses
nanodesign campaign-report campaigns/my-poses
nanodesign campaign-run campaigns/my-poses --max-jobs 1
```

The manifest uses relative paths and remains valid after copying. The calculation stage is a single-point energy/force evaluation on unrelaxed structures. This neither computes a reaction barrier nor selects a physically successful tool. Each job may take many minutes on a CPU; the saved single-geometry timing comparison is under `data/validation`.
