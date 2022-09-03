import pygame
pygame.font.init()

class Button:
	def __init__(self, x, y, width, height, color = (22,59,109), text = None, text_color = (150,255,150), image = None):
		self.x = x
		self.y = y
		self.width = width
		self.height = height
		self.color = color
		self.text = text
		self.text_color = text_color
		self.image = image
		if image != None:
			self.img = pygame.image.load(image).convert_alpha()
			width_ratio = self.width / self.img.get_size()[0]
			height_ratio = self.height / self.img.get_size()[1]
			scale = min(width_ratio, height_ratio)
			img_width = int(self.img.get_size()[0] * scale)
			img_height = int(self.img.get_size()[1] * scale)
			self.img = pygame.transform.scale(self.img, (img_width, img_height))

	def draw(self, screen, size):
		rect = [self.x, self.y, self.width, self.height]
		pygame.draw.rect(screen, self.color, rect)
		if self.text != None:
			font = pygame.font.SysFont("comicsansms", size)
			text_render = font.render(self.text, True, self.text_color)
			text_rect = text_render.get_rect(center = (self.x + self.width/2, self.y + self.height/2))
			screen.blit(text_render, text_rect)
		if self.image != None:
			center_x = self.x + self.width/2
			center_y = self.y + self.height/2
			left = center_x - self.img.get_size()[0]/2
			top = center_y - self.img.get_size()[1]/2
			image_rect = [left, top, self.img.get_size()[0], self.img.get_size()[1]]
			screen.blit(self.img, image_rect)
	
	def contains(self, mouse):
		return self.x < mouse[0] and mouse[0] <= self.x + self.width and self.y < mouse[1] and mouse[1] <= self.y + self.height
