import yt

# Load the dataset.
ds = yt.load("/data/lizli/qso/test/QSO_hdf5_plt_cnt_0130")

print(dir(ds.fields))

for i in sorted(ds.field_list):
  print(i)

for field in ds.fields.gas:
    print(field)

# Create a sphere of radius 100 kpc in the center of the domain.
my_sphere = ds.sphere("c", (500.0, "kpc"))

# Create a PhasePlot object.
# Setting weight to None will calculate a sum.
# Setting weight to a field will calculate an average
# weighted by that field.
plot = yt.PhasePlot(
    my_sphere,
    ("gas", "temperature"),
    ("gas", "density"),
    ("gas", "cell_mass"),
    weight_field=None,
)

#plot.annotate_timestamp(corner='lower_left', time_format='t = {time:.2f} {units}', time_unit="Gyr", text_args={'color':'black'})
plot.set_cmap(("gas", "cell_mass"), "Oranges")

# Set the units of mass to be in solar masses (not the default in cgs)
plot.set_unit(("gas", "cell_mass"), "Msun")
plot.set_font({"size": 35})
plot.set_xlim(1, 1e8)
plot.set_ylim(1e-27, 1e-20)
plot.set_zlim(("gas", "cell_mass"), 1e4, 1e10)

plot.save()
