'''
heightmap.py

python program to apply azimuthal equal area projection to a height map
PyGMT or matplotlib.basemap

python program to generate an stl file

input terrain image map
filename: HuffmanPatterson_2100x1050.png
width: 2100
height:1050

output STL file, for use in 3D printer
filename earth_petri.stl
height 11 mm
diameter 86 mm

project data to azimuthal equal area projection 
center at 1 degree S and 33 degrees E

---------

greek
	Θ θ theta
	Λ λ lambda
	Φ φ phi
	Ψ ψ psi
	Π π pi
	Τ τ tau
	Ρ ρ rho
	°   degrees

angular coordinates can be in degrees or radians, 2π radians = 360 degrees

2D polar coordinates: r, θ
r = radius, distance from origin to point 
θ = angle, angular distance at origin between reference direction and point

3D cylindrical coordinates: ρ, φ, z
ρ = radius, distance along the polar axis, the ray extending outward from the origin in the reference plane
φ = angle at origin between polar axis and point
z = height, distance along longitudinal or cylindrical axis from reference plane to plane of point

3D spherical coordinates: r, θ, φ
r = radius
φ = angle on the horizontal plane from the x axis
θ = angle on the vertical plane from the z axis

geographic coordinates:
λ = longitude
φ = latitude

'''

from PIL import Image
import numpy as np

input_fname = 'earth_heightmap/HuffmanPatterson_2100x1050.png' 

minx	= 0
maxx	= 2100
minlat	= -180.
maxlat	= 180.

miny	= 0
maxy	= 1050
minlon	= -90.
maxlon	= 90.

lon_center = 33.
lat_center = 1.
lon_center = 0.
lat_center = 0.
radius = 1050
lon_center_rad = np.deg2rad(lon_center)
lat_center_rad = np.deg2rad(lat_center)

def interpolate(x,y):
	lat = (x - minx) / (maxx - minx) * (maxlat - minlat) + minlat
	lon = (y - miny) / (maxy - miny) * (maxlon - minlon) + minlon
	return lat,lon

def printHistogram(img_array):
	hist, bin_edges = np.histogram(img_array, bins=256, range=(0, 256))
	print("Frequency Distribution:")
	for i in range(256):
		print(f"Pixel value {i}: {hist[i]}")

def laea_equatorial(lat,lon,prnt=False):  # lambert azimuthal equal area, equatorial
	# equatorial, oblique, polar north, polar south
	# https://github.com/OSGeo/PROJ/blob/master/src/projections/laea.cpp#L102
	# no central point ??

	x = 0.
	y = 0.

	lam = np.deg2rad(lon)
	phi = np.deg2rad(lat)

	sinphi = np.sin(phi)
	sinlam = np.sin(lam)
	cosphi = np.cos(phi)
	coslam = np.cos(lam)

	y = 1. + cosphi * coslam

	if y <= .0000000001:
		print( "y too small")
		return False, False

	y = np.sqrt(2. / y)
	x = y * cosphi * sinlam
	y *= sinphi

	x = int(750 + (x * radius))
	y = int(750 + (y * radius))
	if prnt:
		print(x,y)
	return x,y

def azimuthal_equal_area_wolfram(lat, lon, prnt=False):
	# https://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html
	# no radius?
	λ = np.deg2rad(lon) # lon in radians
	λ0 = np.deg2rad(lon_center) # central longitude in radians
	φ = np.deg2rad(lat) # lat in radians
	φ1 = np.deg2rad(lat_center) # central latitude in radians

	k = np.sqrt(2/(1 + np.sin(φ1) * np.sin(φ) + np.cos(φ1) * np.cos(φ) * np.cos(λ-λ0))) 

	x = k * np.cos(φ) * np.sin(λ-λ0)
	y = k * ((np.cos(φ1) * np.sin(φ)) - (np.sin(φ1) * np.cos(φ) * np.cos(λ-λ0)))

	x = int(radius + (x * radius))
	y = int(radius + (y * radius))

	if prnt:
		print(x,y)
	return x,y

def orthographic_wolfram(lon, lat, prnt=False):
	λ = np.deg2rad(lon)		# lon in radians
	λ0 = np.deg2rad(lon_center)	# central longitude in radians
	φ = np.deg2rad(lat)		# lat in radians
	φ1 = np.deg2rad(lat_center)	# central latitude in radians

	x = np.cos(φ) * np.sin(λ-λ0)
	y = ((np.cos(φ1) * np.sin(φ)) - (np.sin(φ1) * np.cos(φ) * np.cos(λ-λ0)))

	x = int(radius + ((x) * radius))
	y = int(radius + ((y) * radius))
	if prnt:
		print(x,y)
	return x,y

def main():
	global maxx, maxy, radius
	img = Image.open(input_fname).convert('L')  # open as grayscale, default is RGBA
	#img.show()
	
	# Get the dimensions of the image
	width, height = img.size
	maxx, maxy = img.size
	radius = 300 #height
	#print(width,height)
	
	# Convert the image to a numpy array
	img_array = np.array(img)
	#printHistogram(img_array)
	
	# Calculate the new dimensions for the Lambert azimuthal equal-area projection
	new_width = 1500
	new_height = 1500
	
	# Create a new image with the calculated dimensions
	new_img = Image.new('L', (new_width, new_height))  # L:grayscale
	
	# Iterate over the original image and calculate the corresponding coordinates in the new image
	bigx = 0
	bigy = 0
	for x in range(width-1):
		for y in range(height-1):
			lat, lon = interpolate(x,y)
			#new_x, new_y = azimuthal_equal_area_wolfram(lon, lat)
			new_x, new_y = laea_equatorial(lon, lat, False)
			#new_x, new_y = orthographic_wolfram(lon, lat)
			new_img.putpixel((new_x, new_y), int(img_array[y,x]))
			if new_x > bigx:
				bigx = new_x
			if new_y > bigy:
				bigy = new_y

	print(bigx,bigy)
	# Display the converted image
	output_fname = 'earth_heightmap/HuffmanPatterson_lambert_1500x1500.png' 
	output_fname = 'earth_heightmap/HuffmanPatterson_orthographic_1500x1500.png' 
	output_fname = 'earth_heightmap/HuffmanPatterson_laea_eq_1500x1500.png' 
	output_fname = 'earth_heightmap/HuffmanPatterson_lambert_rev_1500x1500.png' 
	new_img.save(output_fname)
	new_img.show()

if __name__ == '__main__':

	#points = [(40.0, 30.0), (50.0, 40.0), (60.0, 50.0)]
	x1, y1 = azimuthal_equal_area_wolfram(33., 1., True)
	x1, y1 = azimuthal_equal_area_wolfram(33., 10., True)
	x1, y1 = azimuthal_equal_area_wolfram(33., -10., True)
	x1, y1 = azimuthal_equal_area_wolfram(23., 1., True)
	x1, y1 = azimuthal_equal_area_wolfram(43., 1., True)

	main()


