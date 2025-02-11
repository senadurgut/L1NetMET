import numpy as np
import pandas as pd
import awkward as ak
import utils.branches as branches
import uproot

def getArrays(inputFiles, branches, nFiles=1, fname="data.parquet"):

    files = [{file: 'Events'} for file in inputFiles][:nFiles]

    # get the data
    data = ak.concatenate([batch for batch in uproot.iterate(files, filter_name=branches)])
    data = formatBranches(data)
    if fname:
        ak.to_parquet(data, fname)

    return data


def getL1Types(useEmu=False, useMP=False):
    
    l1Type = 'L1Emul' if useEmu else 'L1' 
    l1SumType = l1Type + 'MP' if useMP else l1Type
    
    return l1Type, l1SumType


def getBranches(inputs=[], useEmu=False, useMP=False):
    
    l1Type, l1SumType = getL1Types(useEmu, useMP)
    
    sumBranches = [l1SumType + var for var in branches.sumBranches]
    all_branches = sumBranches + branches.puppiMETBranches + branches.muonBranches + branches.recoBranches
    
    for input in inputs:
        all_branches += [l1Type + input + "_" + var for var in branches.objectBranches]

    return all_branches

def formatBranches(data):
    
    # remove the prefixes to the branch names for tidyness
    for branch in ak.fields(data):
        if "L1" in branch:
            data[branch.replace("L1", "").replace("MP", "").replace("Emul", "")] = data[branch]
            del data[branch]
            
    return data
        

def getPUPPIMET(data):
    
    # get the offline puppi MET
    puppiMET = data[branches.puppiMETBranches]
    puppiMET = ak.with_field(puppiMET, puppiMET['PuppiMET_pt']*np.cos(puppiMET['PuppiMET_phi']), "PuppiMET_ptx")
    puppiMET = ak.with_field(puppiMET, puppiMET['PuppiMET_pt']*np.sin(puppiMET['PuppiMET_phi']), "PuppiMET_pty")

    # get the offline muons
    muons = data[branches.muonBranches]
    muons = muons[muons["Muon_isPFcand"] == 1]
    del muons["Muon_isPFcand"]
    muons = ak.with_field(muons, muons['Muon_pt']*np.cos(muons['Muon_phi']), "Muon_ptx")
    muons = ak.with_field(muons, muons['Muon_pt']*np.sin(muons['Muon_phi']), "Muon_pty")

    # make the offline puppi MET no mu
    puppiMET_noMu = ak.copy(puppiMET)
    puppiMET_noMu['PuppiMET_ptx'] = puppiMET['PuppiMET_ptx'] + np.sum(muons['Muon_ptx'], axis=1)
    puppiMET_noMu['PuppiMET_pty'] = puppiMET['PuppiMET_pty'] + np.sum(muons['Muon_pty'], axis=1)
    puppiMET_noMu['PuppiMET_pt'] = np.sqrt(puppiMET_noMu['PuppiMET_ptx']**2 + puppiMET_noMu['PuppiMET_pty']**2)
    
    del puppiMET['PuppiMET_phi'], puppiMET['PuppiMET_ptx'], puppiMET['PuppiMET_pty']
    del puppiMET_noMu['PuppiMET_phi'], puppiMET_noMu['PuppiMET_ptx'], puppiMET_noMu['PuppiMET_pty']
    
    return puppiMET, puppiMET_noMu

def apply_pt_cut(data, puppiMET_noMu, cut_value = -1):
    return data[puppiMET_noMu['PuppiMET_pt'] > cut_value], puppiMET_noMu[puppiMET_noMu['PuppiMET_pt'] > cut_value]

def remove_saturated(data, puppiMET_noMu):
    for col in data.columns:
        if "pt" in col:
            data = data[data[col] < 1000]
            puppiMET_noMu = puppiMET_noMu[data[col] < 1000]
    return data, puppiMET_noMu

def flatten(data, puppiMET_noMu, types=[]):
    
    cutoff = 800
    a = 0.88
    b = 0.06
    c = 0.01

    if 'puppi' in types:
        rand_arr = np.random.rand(len(puppiMET_noMu))
        data = data[rand_arr[puppiMET_noMu['PuppiMET_pt'] > 0]*(a-puppiMET_noMu['PuppiMET_pt']**b/cutoff**b) < c]
        puppiMET_noMu = puppiMET_noMu[rand_arr[puppiMET_noMu['PuppiMET_pt'] > 0]*(a-puppiMET_noMu['PuppiMET_pt']**b/cutoff**b) < c]
    
    if 'l1' in types:
        l1MET = ak.flatten(getSum(data, 'methf')['EtSum_pt'])
        rand_arr = np.random.rand(len(l1MET))
        data = data[rand_arr[l1MET > 0]*(a-l1MET**b/cutoff**b) < c]
        puppiMET_noMu = puppiMET_noMu[rand_arr[l1MET > 0]*(a-l1MET**b/cutoff**b) < c]

    return data, puppiMET_noMu
    

def getCollections(data, inputSums, inputs=[]):

    collections = {}

    # get the sums
    l1Sums = data[branches.sumBranches]
    l1Sums = l1Sums[l1Sums['EtSum_bx'] == 0]
    del l1Sums['EtSum_bx']

    # make the sum collections
    for esum in inputSums:
        sumCol = l1Sums[l1Sums['EtSum_etSumType'] == branches.sums[esum]]
        del sumCol['EtSum_etSumType']
        collections[esum] = sumCol
        
    # make the object collecions
    for input in inputs:
        collection = data[[input + "_" + var for var in branches.objectBranches]]
        collection = collection[collection[input+"_bx"] == 0]
        del collection[input + '_bx']
        collections[input] = collection
        
    return collections

def getSum(data, sumType):
    
    # get the sums
    l1Sums = data[branches.sumBranches]
    l1Sums = l1Sums[l1Sums['EtSum_bx'] == 0]
    del l1Sums['EtSum_bx']

    etSum = l1Sums[l1Sums['EtSum_etSumType'] == branches.sums[sumType]]
    del etSum['EtSum_etSumType']
    
    return etSum


def makeDataframe(collections, fileName=None, nObj=0, keepStruct=False):
    
    object_dfs = []
    for coll in collections:
        if coll in ['Jet', 'EG', 'Tau']:
            objects = pd.DataFrame(ak.to_list(ak.fill_none(ak.pad_none(ak.sort(collections[coll], ascending=False), nObj, clip=True), 0)))
        else:
            objects = pd.DataFrame(ak.to_list(collections[coll]))
        object_labels= ["{}_{}".format(coll, i) for i in range(len(objects.values.tolist()[0][0]))]
        for column in objects.columns:
            object = pd.DataFrame(objects.pop(column).values.tolist())
            object.columns = pd.MultiIndex.from_product([object_labels, [column.split("_")[1]]])
            object_dfs.append(object)

    df = pd.concat(object_dfs, axis=1)

    new_cols = pd.MultiIndex.from_tuples(sorted(list(df.columns)))
    for col in new_cols:
        df[col] = df.pop(col)

    if keepStruct:
        df
    else:
        df.columns = ["{}_{}".format(col[0], col[1]) for col in df.columns]

    if fileName:
        df.to_hdf(fileName, 'online', mode='w')

    return df


def arrayToDataframe(array, label, fileName):

    df = pd.DataFrame(ak.to_list(array))
    if fileName:
        df.to_hdf(fileName, label, mode='a')
    
    return df

def compute_mjjjj(data, puppiMET_noMu):
    def compute_four_jet_mass(pt1, eta1, phi1, pt2, eta2, phi2, pt3, eta3, phi3, pt4, eta4, phi4):
        # Compute individual jet momenta
        def jet_four_momentum(pt, eta, phi):
            px = pt * np.cos(phi)
            py = pt * np.sin(phi)
            pz = pt * np.sinh(eta)
            E  = np.sqrt(px**2 + py**2 + pz**2)  # Assuming massless jets
            return px, py, pz, E

        # Get four-momentum components for each jet
        px1, py1, pz1, E1 = jet_four_momentum(pt1, eta1, phi1)
        px2, py2, pz2, E2 = jet_four_momentum(pt2, eta2, phi2)
        px3, py3, pz3, E3 = jet_four_momentum(pt3, eta3, phi3)
        px4, py4, pz4, E4 = jet_four_momentum(pt4, eta4, phi4)

        # Compute total four-momentum
        E_total  = E1 + E2 + E3 + E4
        px_total = px1 + px2 + px3 + px4
        py_total = py1 + py2 + py3 + py4
        pz_total = pz1 + pz2 + pz3 + pz4

        # Compute invariant mass squared 
        M_inv_squared = np.maximum(E_total**2 - px_total**2 - py_total**2 - pz_total**2, 0)
        return M_inv_squared

    # Compute invariant mass for the four leading jets
    data['m_jjjj_sq'] = compute_four_jet_mass(
        data['Jet_0_pt'], data['Jet_0_eta'], data['Jet_0_phi'],
        data['Jet_1_pt'], data['Jet_1_eta'], data['Jet_1_phi'],
        data['Jet_2_pt'], data['Jet_2_eta'], data['Jet_2_phi'],
        data['Jet_3_pt'], data['Jet_3_eta'], data['Jet_3_phi']
    )
    
    return data, puppiMET_noMu

