from iso_surface import iso_surface

iso = iso_surface() #create_cmap=True)

iso.open('../examples/beltrami_32_fields.nc')

#iso.render(step=60, niso=10)
#iso.export(file_path='./', file_name='test.png')

iso.save_camera_orbiting_animation(step=60, n_frames=360, file_path='./',
                                   file_name="test.mp4", keep_frames=False)


#iso.save_animation(beg=0, end=100,
                   #fps=25,
                   #file_path='./',
                   #file_name="animation.mp4",
                   #keep_frames=True)

iso.close()
