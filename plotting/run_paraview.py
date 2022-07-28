from iso_surface import iso_surface

iso = iso_surface() #create_cmap=True)

iso.open('../examples/beltrami_32_fields.nc')

iso.render(step=60, niso=10)
iso.export(file_path='./', file_name='test.png')


iso.close()
