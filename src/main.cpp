#include "Fluid.h"
#include <SFML/Graphics.hpp>

int main()
{
	//fluid parameters
	int xSize = 200;
	int ySize = 200;
	int elementSize = 4;
	float diff = 0.07;
	float visc = 0.9975;

	//simulation parameters
	int framerate = 60;
	float dt = 1.0 / framerate;
	int numIterations = 5;

	//mouse parameters
	int mouseX = 0;
	int mouseY = 0;
	int lastMouseX = 0;
	int lastMouseY = 0;

	//setting up render window
	sf::RenderWindow window(sf::VideoMode((xSize+2)*elementSize, (ySize+2)*elementSize), "Project");
	window.setFramerateLimit(framerate);

	//creating rendering tools
	sf::Image image;
	sf::Texture texture;
	sf::Sprite sprite;
	image.create((xSize+2)*elementSize, (ySize+2)*elementSize, sf::Color::Black);

	//now let's create the fluid we want (x, y, diffusion, viscosity, timestep, accuracy/quality)
	Fluid fluid(xSize, ySize, diff, visc, dt, numIterations);

	while (window.isOpen())
	{
		//window details
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
				window.close();
		}

		lastMouseX = mouseX;
		lastMouseY = mouseY;
		sf::Vector2i localPosition = sf::Mouse::getPosition(window);
		mouseX = localPosition.x;
		mouseY = localPosition.y;

		if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
			if ((mouseX / elementSize  > 0 && mouseX / elementSize  < xSize - 1 && mouseY / elementSize  > 0 && mouseY / elementSize < ySize - 1) && (lastMouseX / elementSize  > 0 && lastMouseX / elementSize  < xSize - 1 && lastMouseY / elementSize  > 0 && lastMouseY / elementSize < ySize - 1)) {
				fluid.addSource(lastMouseX / elementSize, lastMouseY / elementSize, mouseX / elementSize, mouseY / elementSize);
			}
		}

		if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
			fluid.reset();
		}
		
		fluid.update();
		fluid.render(elementSize, image);
		texture.loadFromImage(image);
		sprite.setTexture(texture, true);
		window.draw(sprite);
		window.display();
	}
    return 0;
}
