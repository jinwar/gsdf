~/bin/gsdf_fit:*.f
	g77 -g -m32 -o ~/bin/gsdf_fit envelop.f fit.f coolb.f inv.f main.f maxsp.f `sac-config --cflags` `sac-config --libs sac` `sac-config --libs sacio`
