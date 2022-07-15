import os, time, random
from turtle import width
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide" #Hide the pygame welcome printing
DIRPATH = os.path.dirname(os.path.abspath(__file__))
os.chdir(DIRPATH) #Sets the directory containing file as the work directory
import pygame, SortingAlgorithms

#Initial Variables
WIDTH = 1000; HEIGHT = 750

MAX = 1000.0; MIN = 0.0; LEN = 200
ALGORITHM = SortingAlgorithms.ShellSort 

CAPTION = "Sorting Algorithms Visualizer - (" + ALGORITHM.__name__ + ")"
ICON = pygame.image.load("ascending.ico")
WIN = pygame.display.set_mode((WIDTH,HEIGHT))

#Display Settings
pygame.display.set_caption(CAPTION)
pygame.display.set_icon(ICON)
pygame.display.flip()


#Functions and Classes
class Bar:
    
    def __init__(self, element, value, win_width, win_height, total_elements, maximum):
        width = win_width/total_elements
        height = win_height/(1.1*maximum)
        self.win_height = win_height
        self.width = width
        self.x = width*element
        self.y = value*height
        self.color = (0, 55 + 200 * value/maximum, 55 + 200 *(1-value/maximum))
        return

    def draw(self, window):
        """Draw the bar"""
        pygame.draw.rect(window, self.color, (self.x, self.win_height - self.y, self.width, self.win_height))
        return



def random_array(max, min, length):
    """Creates a random array with given length, maximum and minimum elements"""
    array = [int((max-min)*random.random()+min) for i in range(length)]
    return array

def draw_array(array, window, max):
    """Draw all the elements in the array in form of bars"""
    window.fill((255,255,255))
    win_width, win_height = window.get_size() #pygame.display.get_surface().get_size()
    for element, value in enumerate(array):
        bar = Bar(element, value, win_width, win_height , len(array), max)
        bar.draw(window)
    pygame.display.update()
    return

def main(window, sorting_method, max, min, length):
    #Create array
    array = random_array(max, min, length)
    draw_array(array, window, max)

    run = True
    while run:
        
        for event in pygame.event.get():

            if event.type == pygame.QUIT: #Quitting the program
                run = False

            if event.type == pygame.KEYDOWN:

                if event.key == pygame.K_SPACE: #If we press space -> Start the sorting
                    start = time.time()
                    copy_array = array.copy()
                    sorting_method(copy_array)
                    run_time = time.time() - start
                    sorting_method(array, lambda: draw_array(array, window, max))
                    vis_time = time.time() - start - run_time
                    print(f"\nSort Algorithm: {sorting_method.__name__}\nSorting Time: {run_time} s.\nVisualization Time: {vis_time} s.")

                if event.key == pygame.K_c: #If we press c -> Reset the array
                    array = random_array(max, min, length)
                    draw_array(array, window, max)


    pygame.quit()

    return array

#Code
if __name__ == '__main__':
    main(WIN, ALGORITHM, MAX, MIN, LEN)