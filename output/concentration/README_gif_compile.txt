ffmpeg -f image2 -framerate 100 -i theta_ink_vol.%04d.png -loop -1 -t 3 theta_ink_vol_100fps.gif



ffmpeg -y -i M_1000fps.gif -filter_complex "fps=5,scale=1080:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=128[p];[s1][p]paletteuse=dither=bayer" output.gif