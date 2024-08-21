class filepaths():
    def __init__(self):
        self.W51dir = '/orange/adamginsburg/w51/TaehwaYoo/'
        #---------------------------------------- W51 B6 (old) ----------------------------------------
        self.W51_b6_cont_old = '/orange/adamginsburg/w51/TaehwaYoo/b6contfits/'
        self.w51e2_b6_briggs = self.W51_b6_cont_old+'W51e2_cont_bigbriggs.image.fits'
        self.w51e2_b6_robust0 = self.W51_b6_cont_old+'W51e2_cont_big_robust0.image.fits'
        self.w51e2_b6_uniform = self.W51_b6_cont_old+'W51e2_cont_biguniform.image.fits'
        self.w51e2_b6_superuniform = self.W51_b6_cont_old+'W51e2_cont_bigsuperuniform.image.fits'

        self.w51n_b6_briggs = self.W51_b6_cont_old+'W51n_cont_bigbriggs.image.fits'
        self.w51n_b6_robust0 = self.W51_b6_cont_old+'w51n_cont_big_robust0.image.fits'
        self.w51n_b6_uniform = self.W51_b6_cont_old+'W51n_cont_biguniform.image.fits'
        self.w51n_b6_superuniform = self.W51_b6_cont_old+'W51n_cont_bigsuperuniform.image.fits'
        self.w51n_b6_natural = self.W51_b6_cont_old+'W51n_cont_bignatural.image.fits'

        #---------------------------------------- W51 B6 (new) ----------------------------------------
        self.w51e_b6_cont = '/orange/adamginsburg/w51/TaehwaYoo/w51e2.spw0thru19.14500.robust0.thr0.15mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
        self.w51n_b6_cont = '/orange/adamginsburg/w51/TaehwaYoo/w51n.spw0thru19.14500.robust0.thr0.1mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'

        #---------------------------------------- W51 B3 ----------------------------------------
        self.W51b3 = '/orange/adamginsburg/w51/TaehwaYoo/2017.1.00293.S_W51_B3_LB/may2021_successful_imaging/'
        self.w51n_b3_tt0 = self.W51b3+'w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
        self.w51n_b3_tt1 = self.W51dir +'w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt1.pbcor.fits'
        self.w51n_b3_alpha = self.W51dir +'w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.alpha.pbcor.fits'


        self.w51e_b3_tt0 = self.W51b3+'w51e2.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
        self.w51e_b3_tt1 = self.W51dir+'w51e2.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt1.pbcor.fits'
        self.w51e_b3_alpha = self.W51dir+'w51e2.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.alpha.pbcor.fits'


        #---------------------------------------- image convolved to the common beam size ----------------------------------------
        self.w51conv = '/orange/adamginsburg/w51/TaehwaYoo/convolved_new/'
        self.w51e_b6_conv = self.w51conv + 'w51e_new_B6_conv.fits'
        self.w51n_b6_conv = self.w51conv + 'w51n_new_B6_conv.fits'
      
        self.w51e_b6_almaimf_conv = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/w51e_B6_conv.fits'
        self.w51e_b3_almaimf_conv = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/w51e_B3_conv.fits'
        self.w51n_b3_almaimf_conv = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/w51n_B3_conv.fits'
        self.w51n_b6_almaimf_conv = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/w51n_B6_conv.fits'

        
        #---------------------------------------- ALMA-IMF W51 ----------------------------------------

        self.w51n_b6_almaimf = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/W51-IRS2_B6_uid___A001_X1296_X187_continuum_merged_12M_robust0_selfcal9_finaliter.image.tt0.pbcor.fits'

        #---------------------------------------- W51 background noises ----------------------------------------
        self.w51e_b6_background = '/home/t.yoo/w51/catalogue/dendrogram/background_noise_w51e_b6_pb.dat'
        self.w51e_b3_background = '/home/t.yoo/w51/catalogue/dendrogram/background_noise_w51e_b3_pb.dat'
        self.w51n_b6_background = '/home/t.yoo/w51/catalogue/dendrogram/background_noise_w51n_b6_pb.dat'
        self.w51n_b3_background = '/home/t.yoo/w51/catalogue/dendrogram/background_noise_w51n_b3_pb.dat'


        #---------------------------------------- W51 catalogs ----------------------------------------
        self.w51e_dendro_matched_catalog = '/home/t.yoo/w51/catalogue/dendrogram/dendro_w51e_matched.fits'
        self.w51n_dendro_matched_catalog = '/home/t.yoo/w51/catalogue/dendrogram/dendro_w51n_matched.fits'

        #---------------------------------------- localdir ----------------------------------------
        self.localdir = '/Users/dbahck37/w51data'
        self.w51e_b6_cont_local = self.localdir + '/w51e2.spw0thru19.14500.robust0.thr0.15mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
        self.w51n_b6_cont_local = self.localdir + '/w51n.spw0thru19.14500.robust0.thr0.1mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
        self.w51e_b3_cont_local = self.localdir + '/b3data/w51e2.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
        self.w51n_b3_cont_local = self.localdir + '/b3data/w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
        self.w51e_b6_conv_local = self.localdir + '/convolved_new/w51e_new_B6_conv.fits'
        self.w51n_b6_conv_local = self.localdir + '/convolved_new/w51n_new_B6_conv.fits'    

        self.w51e_dendro_matched_catalog_local = '/Users/dbahck37/w51_jupyter/w51/catalogue/dendrogram/dendro_w51e_matched.fits'
        self.w51n_dendro_matched_catalog_local = '/Users/dbahck37/w51_jupyter/w51/catalogue/dendrogram/dendro_w51n_matched.fits'

