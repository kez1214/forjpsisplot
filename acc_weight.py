#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT
import glob
from array import array
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ==================================================
# Input files
# ==================================================
files = glob.glob(
    "/eos/user/z/zhipeng/MC_Particle_Gun_For_Acc/results/step0_*.root"
)
if len(files) == 0:
    raise RuntimeError("No input files found!")
print("Total files =", len(files))

refFile = "fitOutput.root"

# ==================================================
# Read reference distribution from fitOutput.root
# ==================================================
fRef = ROOT.TFile.Open(refFile)
hRef = fRef.Get("hPromptPtY")
if not hRef:
    raise RuntimeError("Cannot find hPromptPtY")

print(f">>> hRef entries = {hRef.GetEntries()}, integral = {hRef.Integral():.2f}")

# ==================================================
# Binning
# ==================================================
y_bins = array('d', [0.0, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4])
pt_bins = array('d', [
    10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    32, 34, 36, 38, 42, 46, 50, 60, 75, 95, 120, 200
])

nY = len(y_bins) - 1
nP = len(pt_bins) - 1

# ==================================================
# Histograms
# ==================================================
hGenPtY = ROOT.TH2D("hGenPtY", "Gen J/psi;p_{T};|y|", nP, pt_bins, nY, y_bins)

hDenom = ROOT.TH2D("hDenom", "Denominator;|y|;p_{T} [GeV]", nY, y_bins, nP, pt_bins)
hNumer = ROOT.TH2D("hNumer", "Numerator;|y|;p_{T} [GeV]",   nY, y_bins, nP, pt_bins)

hDenomW = ROOT.TH2D("hDenomW", "Weighted Denominator;|y|;p_{T} [GeV]", nY, y_bins, nP, pt_bins)
hNumerW = ROOT.TH2D("hNumerW", "Weighted Numerator;|y|;p_{T} [GeV]",   nY, y_bins, nP, pt_bins)

# ==================================================
# Acceptance criteria
# ==================================================
def muon_pass_criteria(mu):
    return abs(mu.eta()) < 2.4 and mu.pt() > 5.0

# ==================================================
# FWLite setup
# ==================================================
handle = Handle("std::vector<reco::GenParticle>")
label  = "genParticles"

# ==================================================
# FIRST PASS: fill gen pT-|y| distribution
# ==================================================
print(">>> First pass: filling gen distribution")
events = Events(files)

nEvents = 0
nJpsi_gen = 0

for iev, event in enumerate(events):
    if iev % 10000 == 0:
        print("  first pass event =", iev)
    nEvents += 1

    event.getByLabel(label, handle)
    gen_particles = handle.product()

    for p in gen_particles:
        if p.pdgId() != 443:
            continue

        if not p.isLastCopy():
            continue

        pt = p.pt()
        y  = abs(p.rapidity())

        # 只保留在 bin 范围内的事件
        if not (pt_bins[0] <= pt < pt_bins[-1]):
            continue
        if not (y_bins[0] <= y < y_bins[-1]):
            continue

        hGenPtY.Fill(pt, y)
        nJpsi_gen += 1

print(f">>> First pass: events = {nEvents}, gen J/psi in range = {nJpsi_gen}")

# ==================================================
# Compute reweighting map
# ==================================================
hWeight = hRef.Clone("hWeight")
hWeight.Reset()
hWeight.SetTitle("Reweighting map;p_{T};|y|")

refInt = hRef.Integral()
genInt = hGenPtY.Integral()
print(f">>> refInt = {refInt:.2f}, genInt = {genInt:.2f}")

for ix in range(1, hRef.GetNbinsX() + 1):
    for iy in range(1, hRef.GetNbinsY() + 1):
        refFrac = hRef.GetBinContent(ix, iy) / refInt if refInt > 0 else 0.0
        genFrac = hGenPtY.GetBinContent(ix, iy) / genInt if genInt > 0 else 0.0

        if genFrac > 0:
            hWeight.SetBinContent(ix, iy, refFrac / genFrac)
        else:
            hWeight.SetBinContent(ix, iy, 0.0)

# ==================================================
# Weight function
# ==================================================
def getWeight(pt, y):
    y = abs(y)
    if not (pt_bins[0] <= pt < pt_bins[-1]):
        return 0.0
    if not (y_bins[0] <= y < y_bins[-1]):
        return 0.0

    ix = hWeight.GetXaxis().FindFixBin(pt)
    iy = hWeight.GetYaxis().FindFixBin(y)
    return hWeight.GetBinContent(ix, iy)

# ==================================================
# SECOND PASS: fill acceptance histograms with and without weights
# ==================================================
print(">>> Second pass: filling acceptance histograms")
events = Events(files)

nPassed = 0
sumW = 0.0

for iev, event in enumerate(events):
    if iev % 10000 == 0:
        print("  second pass event =", iev)

    event.getByLabel(label, handle)
    gen_particles = handle.product()

    for p in gen_particles:
        if p.pdgId() != 443:
            continue
        if not p.statusFlags().isLastCopy():
            continue
#lastcopy

        pt = p.pt()
        y  = abs(p.rapidity())

        if not (pt_bins[0] <= pt < pt_bins[-1]):
            continue
        if not (y_bins[0] <= y < y_bins[-1]):
            continue

        w = getWeight(pt, y)

        # 找 muon 子粒子
        daughters = []
        for iD in range(p.numberOfDaughters()):
            dau = p.daughter(iD)
            if abs(dau.pdgId()) == 13:
                daughters.append(dau)

        passed = (len(daughters) == 2 and
                  all(muon_pass_criteria(mu) for mu in daughters))

        # 不加权
        hDenom.Fill(y, pt)
        if passed:
            hNumer.Fill(y, pt)

        # 加权
        hDenomW.Fill(y, pt, w)
        if passed:
            hNumerW.Fill(y, pt, w)

        if passed:
            nPassed += 1
        sumW += w

print(f">>> Second pass: passed = {nPassed}, sum of weights = {sumW:.2f}")

# ==================================================
# Compute acceptance (unweighted and weighted)
# ==================================================
hEff = hNumer.Clone("hEff")
hEff.Divide(hDenom)

hEffW = hNumerW.Clone("hEffW")
hEffW.Divide(hDenomW)

# 设置误差
def set_bin_errors(hEff, hNum, hDen):
    for ix in range(1, hEff.GetNbinsX() + 1):
        for iy in range(1, hEff.GetNbinsY() + 1):
            n_num = hNum.GetBinContent(ix, iy)
            n_den = hDen.GetBinContent(ix, iy)
            if n_den > 0 and 0 <= n_num <= n_den:
                acc = n_num / n_den
                err = (acc * (1.0 - acc) / n_den) ** 0.5
                hEff.SetBinError(ix, iy, err)
            else:
                hEff.SetBinError(ix, iy, 0.0)

set_bin_errors(hEff,  hNumer,  hDenom)
set_bin_errors(hEffW, hNumerW, hDenomW)

# ==================================================
# Save ROOT file
# ==================================================
outFile = ROOT.TFile("output_Acceptance.root", "RECREATE")
hGenPtY.Write()
hWeight.Write()
hDenom.Write()
hNumer.Write()
hEff.Write()
hDenomW.Write()
hNumerW.Write()
hEffW.Write()
outFile.Close()

print(">>> output_Acceptance.root written")

# ==================================================
# Save accData_cpp.txt
# ==================================================
with open("accData_cpp.txt", "w") as out:
    nbins_x = hEffW.GetXaxis().GetNbins()   # |y|
    nbins_y = hEffW.GetYaxis().GetNbins()   # pT

    for iy in range(1, nbins_y + 1):
        pt_low  = hEffW.GetYaxis().GetBinLowEdge(iy)
        pt_high = hEffW.GetYaxis().GetBinUpEdge(iy)
        out.write(f"// pT {pt_low:g}-{pt_high:g}\n")

        for ix in range(1, nbins_x + 1):
            y_low  = hEffW.GetXaxis().GetBinLowEdge(ix)
            y_high = hEffW.GetXaxis().GetBinUpEdge(ix)

            acc     = hEffW.GetBinContent(ix, iy)
            acc_err = hEffW.GetBinError(ix, iy)

            out.write(
                f"{{{y_low:.1f}, {y_high:.1f}, "
                f"{pt_low:.1f}, {pt_high:.1f}, "
                f"{acc:.6f}, {acc_err:.6f}}},\n"
            )
        out.write("\n")

print(">>> accData_cpp.txt saved")

# ==================================================
# Save numerator / denominator (before and after weighting)
# ==================================================
nFiles = len(files)
nTotEvents   = nEvents
nTotGen      = nJpsi_gen
nSelected    = int(hDenom.Integral())
nAccepted    = int(hNumer.Integral())

sumDenW  = hDenomW.Integral()
sumNumW  = hNumerW.Integral()
sumAllW  = sumDenW
intAcc   = (sumNumW / sumDenW) if sumDenW > 0 else 0.0
intAccUn = (nAccepted / nSelected) if nSelected > 0 else 0.0

with open("accNumerDenom.txt", "w") as out:
    out.write("============================================================\n")
    out.write("Numerator / Denominator  (unweighted and weighted)\n")
    out.write("============================================================\n")
    out.write(f"Number of input files       = {nFiles}\n")
    out.write(f"Total events                 = {nTotEvents}\n")
    out.write(f"Total GEN J/psi              = {nTotGen}\n")
    out.write(f"Selected J/psi               = {nSelected}\n")
    out.write(f"Accepted J/psi               = {nAccepted}\n")
    out.write(f"Weighted denominator         = {sumDenW:.10f}\n")
    out.write(f"Weighted numerator           = {sumNumW:.10f}\n")
    out.write(f"Sum of weights               = {sumAllW:.10f}\n")
    out.write(f"Integrated acceptance (unw)  = {intAccUn:.10f}\n")
    out.write(f"Integrated acceptance (w)    = {intAcc:.10f}\n")
    out.write("\n")
    out.write("============================================================\n")
    out.write("Bin-by-bin numerator and denominator\n")
    out.write("============================================================\n")
    out.write(
        f"{'pT_low':<12}{'pT_high':<13}"
        f"{'y_low':<13}{'y_high':<13}"
        f"{'Denominator':<20}{'Numerator':<20}{'Acceptance':<20}"
        f"{'Den_w':<20}{'Num_w':<20}{'Acc_w':<20}\n"
    )
    out.write("-" * 168 + "\n")

    nbins_x = hDenom.GetXaxis().GetNbins()   # |y|
    nbins_y = hDenom.GetYaxis().GetNbins()   # pT

    for iy in range(1, nbins_y + 1):
        pt_low  = hDenom.GetYaxis().GetBinLowEdge(iy)
        pt_high = hDenom.GetYaxis().GetBinUpEdge(iy)

        for ix in range(1, nbins_x + 1):
            y_low  = hDenom.GetXaxis().GetBinLowEdge(ix)
            y_high = hDenom.GetXaxis().GetBinUpEdge(ix)

            den   = hDenom.GetBinContent(ix, iy)
            num   = hNumer.GetBinContent(ix, iy)
            acc   = hEff.GetBinContent(ix, iy)

            den_w = hDenomW.GetBinContent(ix, iy)
            num_w = hNumerW.GetBinContent(ix, iy)
            acc_w = hEffW.GetBinContent(ix, iy)

            out.write(
                f"{pt_low:<12.3f}{pt_high:<13.3f}"
                f"{y_low:<13.3f}{y_high:<13.3f}"
                f"{den:<20.8f}{num:<20.8f}{acc:<20.8f}"
                f"{den_w:<20.8f}{num_w:<20.8f}{acc_w:<20.8f}\n"
            )

    out.write("\n")
    out.write("============================================================\n")
    out.write("Acceptance range check\n")
    out.write("============================================================\n")

    ok = True
    for ix in range(1, nbins_x + 1):
        for iy in range(1, nbins_y + 1):
            a  = hEff.GetBinContent(ix, iy)
            aw = hEffW.GetBinContent(ix, iy)
            if not (0.0 <= a <= 1.0) or not (0.0 <= aw <= 1.0):
                ok = False

    if ok:
        out.write("All acceptance bins are within [0, 1]\n")
    else:
        out.write("WARNING: Some acceptance bins are outside [0, 1]\n")

print(">>> accNumerDenom.txt saved")

# ==================================================
# Draw acceptance maps
# ==================================================
ROOT.gStyle.SetOptStat(0)

# 不加权 acceptance
c1 = ROOT.TCanvas("c1", "Acceptance unweighted", 900, 700)
hEff.SetTitle("Acceptance (unweighted);|y|;p_{T} [GeV]")
hEff.SetMarkerSize(0.8)
hEff.Draw("COLZ TEXT")
c1.SaveAs("AcceptanceMap.png")
print(">>> AcceptanceMap.png saved")

# 加权 acceptance
c2 = ROOT.TCanvas("c2", "Acceptance weighted", 900, 700)
hEffW.SetTitle("Acceptance (weighted by prompt J/#psi data);|y|;p_{T} [GeV]")
hEffW.SetMarkerSize(0.8)
hEffW.Draw("COLZ TEXT")
c2.SaveAs("AcceptanceMap_weighted.png")
print(">>> AcceptanceMap_weighted.png saved")

# 权重图
c3 = ROOT.TCanvas("c3", "Weight map", 900, 700)
hWeight.SetTitle("Reweighting map;p_{T};|y|")
hWeight.SetMarkerSize(0.8)
hWeight.Draw("COLZ TEXT")
c3.SaveAs("WeightMap.png")
print(">>> WeightMap.png saved")

print("===================================")
print("All done.")
print("===================================")
