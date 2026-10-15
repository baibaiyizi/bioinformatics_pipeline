"""科学绘图输入边界：缺测、伪造效应、身份及重复键。"""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import mechanism_recipes as plot

class MechanismContract(unittest.TestCase):
    def cpg(self):
        return pd.DataFrame({'locus_id':['a','a'],'context':['sperm','liver'],
            'sample_id':['s1','l1'],'condition':['C','T'],'chr':['chr1','chr1'],
            'position':[100,100],'methylated_reads':[0,np.nan],'total_reads':[0,np.nan]})

    def test_zero_and_censored_reads_remain_missing(self):
        fig=plot.cpg(self.cpg())
        values=fig.axes[0].images[0].get_array()
        self.assertTrue(np.ma.getmaskarray(values).all())
        plot.plt.close(fig)

    def test_fractional_or_excess_methylated_reads_rejected(self):
        for m,n in [(3,2),(.5,2),(-1,2)]:
            d=self.cpg();d.loc[0,['methylated_reads','total_reads']]=[m,n]
            with self.assertRaises(ValueError):plot.cpg(d)

    def test_duplicate_cpg_library_keys_rejected(self):
        d=self.cpg()
        with self.assertRaises(ValueError):plot.cpg(pd.concat([d,d.iloc[:1]]))

    def test_ineligible_methylation_effect_rejected(self):
        d=pd.DataFrame({'edge_id':['a'],'tf':['X'],'target':['Y'],'motif_label':['M'],
            'methylation_qualified':[False],'tf_activity_delta':[1],
            'target_rna_log2fc':[.2],'motif_deviation_delta':[.3],'methylation_delta':[.4]})
        with self.assertRaises(ValueError):plot.regulatory(d)

    def test_unmeasured_metabolite_effect_rejected(self):
        d=pd.DataFrame({'node_id':['a','b'],'label':['A','B'],'x':[0,1],'y':[0,0],
            'measurement_state':['unmeasured','unmeasured'],'effect':[1,np.nan]})
        e=pd.DataFrame({'source':['a'],'target':['b'],'label':['E'],'edge_type':['reaction'],'source_url':['https://www.rhea-db.org/rhea/64720']})
        with self.assertRaises(ValueError):plot.reactions(d,e)

    def test_unsourced_reaction_rejected(self):
        d=pd.DataFrame({'node_id':['a','b'],'label':['A','B'],'x':[0,1],'y':[0,0],
            'measurement_state':['unmeasured','unmeasured'],'effect':[np.nan,np.nan]})
        e=pd.DataFrame({'source':['a'],'target':['b'],'label':['E'],'edge_type':['reaction'],'source_url':['']})
        with self.assertRaises(ValueError):plot.reactions(d,e)

if __name__=='__main__':unittest.main()
