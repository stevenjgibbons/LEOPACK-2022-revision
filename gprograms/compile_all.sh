#!/bin/sh
for program in \
	arrows_const_r3 \
	arrows_z_eq_merid4 \
	cic_arrows_z_eq_merid4 \
	cicm_arrows_z_eq_merid4 \
	continent_arrows_const_r3 \
	continent_full_sphere_plot \
	cutout_sphere_plot \
	full_sphere_plot \
	ps_2plot_z_eq_merid2 \
	shc_sphere_plot
do
	make $program
done
