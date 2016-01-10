#!/bin/bash

cd ~/galfacts/inpainting/

python ~/repos/galfacts-python/aps/aps_noise_determination.py N2_Q_inpainted.fits N2_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits N2_P_inpainted.fits N2 900 0 2

python ~/repos/galfacts-python/aps/aps_noise_determination.py N3_Q_inpainted.fits N3_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits N3_P_inpainted.fits N3 900 0 2

python ~/repos/galfacts-python/aps/aps_noise_determination.py N4_Q_inpainted.fits N4_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits N4_P_inpainted.fits N4 900 0 2

python ~/repos/galfacts-python/aps/aps_noise_determination.py S1_Q_inpainted.fits S1_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits S1_P_inpainted.fits S1 900 0 2

python ~/repos/galfacts-python/aps/aps_noise_determination.py S2_Q_inpainted.fits S2_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits S2_P_inpainted.fits S2 900 0 2

python ~/repos/galfacts-python/aps/aps_noise_determination.py S3_Q_inpainted.fits S3_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits S3_P_inpainted.fits S3 900 0 2

python ~/repos/galfacts-python/aps/aps_noise_determination.py S4_Q_inpainted.fits S4_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits S4_P_inpainted.fits S4 900 0 2


