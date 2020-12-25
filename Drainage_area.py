from parameters import *

def f_area(area,direction,channel_loc):
	for x in xrange(1,grid_cells+1):
		for y in xrange (2,grid_cells):
			if channel_loc[x][y]==1:
				area_dis=1
			else:
				area_dis=1
			area[x][y]=area[x][y]+1
			xloc = x
			yloc = y
			bingo = 0
			while bingo ==0:
				if yloc ==1 or yloc == grid_cells:
					bingo =1
				else:
					i=direction[xloc][yloc]
					xloc = xloc + xn[i]
					yloc = yloc +yn[i]
					area[xloc][yloc] = area[xloc][yloc]+area_dis
			
	return (area)