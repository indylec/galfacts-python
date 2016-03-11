#!/bin/bash

cd /Users/leclercq/galfacts/inpainting

python ~/repos/galfacts-python/aps/aps_final_results.py N2_Q_inpainted.fits N2_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits 2 N2

python ~/repos/galfacts-python/aps/aps_final_results.py N3_Q_inpainted.fits N3_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits 2 N3

python ~/repos/galfacts-python/aps/aps_final_results.py N4_Q_inpainted.fits N4_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits 2 N4

python ~/repos/galfacts-python/aps/aps_final_results.py S1_Q_inpainted.fits S1_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits 2 S1

python ~/repos/galfacts-python/aps/aps_final_results.py S2_Q_inpainted.fits S2_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits 2 S2

python ~/repos/galfacts-python/aps/aps_final_results.py S3_Q_inpainted.fits S3_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits 2 S3

python ~/repos/galfacts-python/aps/aps_final_results.py S4_Q_inpainted.fits S4_U_inpainted.fits simu_noise_test3_S2_average_image_Q.fits simu_noise_test3_S2_average_image_U.fits 2 S4
