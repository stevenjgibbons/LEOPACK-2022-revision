#
# Makefile for the graphics subroutines
#
PGSOURCES= \
   hls_cont_set.f                 \
   pgconf_polar.f                 \
   polar_sub.f                    \
   const_rad_arc_draw.f           \
   const_the_line_draw.f          \
   cont_fill_polar.f              \
   cont_draw_polar.f              \
   polar_box_draw.f               \
   array_polar_colour.f           \
   array_polar_draw.f             \
   pgplot_file_init.f             \
   extreme_coords_find.f          \
   merid_sec_polar_eval.f         \
   eq_sec_polar_eval.f            \
   constant_z_polar_eval.f        \
   constant_r_rect_eval.f         \
   array_rect_colour.f            \
   array_rect_draw.f              \
   cont_fill_rect.f               \
   cont_draw_rect.f               \
   vs_from_vr_vt.f                \
   array_polar_arrows.f           \
   array_polar_arrows2.f          \
   array_rect_arrows.f            \
   osc_func_eval.f                \
   osc_func_eval_option.f         \
   array_rect_colour_option.f     \
   array_rect_draw_option.f       \
   array_polar_colour_option.f    \
   array_polar_draw_option.f      \
   valmin_valmax_update.f         \
   spher_sat_2_world.f            \
   latdeg_2_the.f                 \
   longdeg_2_phi.f                \
   pgconf_general.f               \
   array_general_colour.f         \
   pgconf_spher_proj.f            \
   array_spher_proj_colour.f      \
   array_spher_proj_colour2.f     \
   array_annihilate.f             \
   spherical_line_draw.f          \
   cut_out_sphere_plot.f          \
   cut_out_sphere_plot2.f         \
   inhomog_temp_vr_temp_eq.f      \
   simple_sphere_plot.f           \
   cont_draw_polar_mod.f          \
   cont_draw_rect_mod.f           \
   array_polar_draw_option_mod.f  \
   array_rect_draw_option_mod.f
#
SOURCES= \
   $(PGSOURCES)
#
OPTIM= -O3
PGSLIB= ../lib/pgslib.a
FORTRAN= gfortran
.PRECIOUS:      $(PGSLIB)
#
$(PGSLIB):       $(SOURCES)
	$(FORTRAN) $(OPTIM) $? -c
	ar rv $(PGSLIB) *.o
	ranlib $(PGSLIB)
	chmod 644 $(PGSLIB)
	\rm *.o
#
