#SETUP: Make sure to install python3. You might have success with regular python, but I haven't tested
#Then install matplotlib. You also need numpy, but I think that is installed by default.
#Then you can run this script with the terminal.
#Navigate to the folder containing this file, and type "python3 braid_visualization.py"

#Here is all the code.
#Scroll to the bottom to see examples and how to draw some braids!

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

class Braid():
    def __init__(self, crossings, n=None):
        self.crossings = self.valid_crossings(crossings)
        self.n = self.valid_number_of_strands(n)
        #we keep track of strand positions to properly color the strands
        self.strand_positions = self.compute_strand_paths()

    def invert(self):
        inverted_crossings = list(map(lambda x: x * -1, self.crossings))[::-1]
        return Braid(inverted_crossings, n=self.n)

    def valid_number_of_strands(self, n):
        largest_sigma = max([abs(i) for i in self.crossings])
        if n == None:
            n = largest_sigma + 1
        elif n < largest_sigma:
            raise ValueError('It is impossible to have {} strands with an S-{} crossing!'.format(n, largest_sigma))
        return n

    def valid_crossings(self, crossings):
        if not all(isinstance(i, int) for i in crossings):
            raise ValueError('All crossings must be integers!\n{}'.format(crossings))
        return crossings

    def draw(self, length = 8, line_width = 4.0, save_name='', show=True):
        plt.figure(figsize=(length, length / 12 * len(self.strand_positions)))
        plt.axis('off')
        colors = [cm.jet(i) for i in np.linspace(0.0, 1.0, len(self.strand_positions))]
        for x in range(len(self.crossings)):
            #separate between crossings and non-crossings
            if self.crossings[x] == 0:
                non_crossing_indices = list(range(self.n))
            else:
                crossing_indices = [abs(self.crossings[x])-1, abs(self.crossings[x])]
                non_crossing_indices = list(range(crossing_indices[0])) + list(range(crossing_indices[1]+1, self.n))
                over  = crossing_indices[0] if np.sign(self.crossings[x]) == -1 else crossing_indices[1]
                under = crossing_indices[0] if np.sign(self.crossings[x]) ==  1 else crossing_indices[1]
                over_data  = self.segment_drawing_data(x, over , -1 * np.sign(self.crossings[x]))
                under_data = self.segment_drawing_data(x, under,      np.sign(self.crossings[x]))
                strand_at_over  = self.pos_of_strand_at(over,  strand_paths=self.strand_positions, x=x)
                strand_at_under = self.pos_of_strand_at(under, strand_paths=self.strand_positions, x=x)
                plt.plot(*over_data , color=colors[strand_at_over] , linewidth=line_width)
                plt.plot(*under_data, color=colors[strand_at_under], linewidth=line_width)
            for i in non_crossing_indices:
                strand_at_i = self.pos_of_strand_at(i, strand_paths=self.strand_positions, x=x)
                data = self.segment_drawing_data(x, i, 0)
                plt.plot(*data, color=colors[strand_at_i], linewidth=line_width)

        if save_name != '':
            plt.savefig(save_name, bbox_inches='tight')
        if show:
            plt.show()

    #given an x (crossing) and y (height on plot), find the index of the strand in this position
    def pos_of_strand_at(self, y, strand_paths, x=-1): #x defaults to most recent strand position
        if not y in range(self.n):
            raise ValueError('It is impossible to have a strand at position {}, when there are only {} positions.'.format(y,self.n))
        #now we don't need to check if there is a strand at the desired position, because there are n strands in n different positions
        #this means all positions must have a strand
        strand_positions_at_x = [i[x] for i in strand_paths]
        return strand_positions_at_x.index(y)


    def compute_strand_paths(self):
        strand_paths = [[i] for i in range(self.n)] #each list contains coordinates of each strand at each time step
        for c in self.crossings:
            #add coordinates for new timestep
            #assume no changes are made and add previous coordinates as the next ones
            for s in strand_paths:
                s.append(s[-1])

            #if there is a crossing, then update the position of the two strands that crossed
            if c != 0:
                lower_strand = self.pos_of_strand_at(abs(c)-1, strand_paths=strand_paths)
                upper_strand = self.pos_of_strand_at(abs(c)  , strand_paths=strand_paths)
                strand_paths[lower_strand][-1] += 1
                strand_paths[upper_strand][-1] -= 1

        return strand_paths

    #x_start should be an integer corresponding to a crossing
    #y_start should be an integer corresponding to the height of a strand at position x
    #direction should be -1, 0, or 1, indicating the end height of y_data
    def segment_drawing_data(self, x_start, y_start, direction, n = 100):
        x_data = np.linspace(x_start, x_start+1, num = n)
        cos_range = np.linspace(0, np.pi, num = n)
        y_data = list(map(lambda x: direction * (1 - np.cos(x)) / 2 + y_start, cos_range))
        return [x_data, y_data]

#this is outside of the class on purpose!
def mult_braids(braids):
    if len(braids) == 0:
        return None
    if any([i.n != braids[0].n for i in braids]):
        raise ValueError('Braids have differing numbers of strands! {}'.format([i.n for i in braids]))
    return Braid(sum([i.crossings for i in braids], []), n=braids[0].n)


'''
Great Jacob, you've made all this code but what can I do with it?
I'm so glad you've asked!
This program only draws braid groups, it doesn't really do anything else.

First, to define a braid, you need to describe it with artin generators!
In most books this is a product of sigma_i's.
For the sake of your poor fingers, instead of sigma_i, just write i for short.
Instead of multiplying these i's, place them in a list!
Instead of writing inverses as i^-1, simply write -i.
So the braid sigma_1, sigma_2^(-1), sigma_3 is written as [1, -2, 3]
You could draw [1, -2, 3] like this.
'''

#uncomment this next line to actually draw it!
#Braid([1, -2, 3]).draw()
#You can save the image with the save_name parameter in draw()
#Just pick the filename you want. For example,
#Braid([[1, -2, 3]).draw(save_name='put_a_file_name_here')
#If you want to save a braid but not show the image, use show=False
#Braid([1, -2, 3]).draw(save_name='this_is_not_shown', show=False)

'''
There are two more important things to note!
One is that instead of an artin generator, you can write 0.
This just adds horizontal space in the image (equal to the space taken by a crossing)
I used this to show visual similarities between different braids.
If you wanted to show that [1, 3, -1] = [3],
you could break it into a few steps that are easier to visualize, like this:
[1, 3, -1] = [1, -1, 3] = [0, 0, 3] = [3]

The final thing to note is that this program draws one braid at a time.
So in the following examples, until you close out of the first drawing, the next doesn't appear.
This gives a flip book like effect which is the next best thing before creating full-blown animation.
'''


#'''
#uncomment to show an example of combing
combing_crossings = []
combing_crossings.append([ 1,-2, 1, 1, 0, 0, 0, 0, 2,-1])
combing_crossings.append([ 1,-2, 1,-2, 0, 0, 0, 0, 1, 2])
combing_crossings.append([ 1,-2, 1,-2,-1, 1, 0, 0, 1, 2])
combing_crossings.append([ 1,-2,-2,-1, 2, 1, 0, 0, 1, 2])
combing_crossings.append([ 1,-2,-2,-1, 2, 1,-2, 2, 1, 2])
combing_crossings.append([ 1,-2,-2,-1,-1, 2, 1, 2, 1, 2])
combing_crossings.append([ 1,-2,-2,-1,-1, 2, 2, 1, 2, 2])
for i in range(len(combing_crossings)):
    Braid(combing_crossings[i]).draw()
#'''


'''
#uncomment to show an example of handle reduction
handle_reduction_crossings.append([1,2,3,-2,-1])
handle_reduction_crossings.append([1,-3,2,3,-1])
handle_reduction_crossings.append([-3,1,2,-1,3])
handle_reduction_crossings.append([-3,-2,1,2,3])
for i in range(len(handle_reduction_crossings)):
    Braid(handle_reduction_crossings[i]).draw(save_name='handle_reduction_{}'.format(i+1))
'''
