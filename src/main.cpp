#include <SFML/Graphics.hpp>
#include <vector>

void drawWaves(sf::RenderWindow& window, const std::vector<std::vector<std::vector<float>>>& amplitudes) {
    // Размер окна
    const int width = window.getSize().x;
    const int height = window.getSize().y;

    // Создаем изображение
    sf::Image image;
    image.create(width, height, sf::Color::White);

    // Проходим по каждой точке изображения и устанавливаем цвет в зависимости от амплитуды
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            float amplitude = amplitudes[x][y][0]; // Здесь используется только одна компонента амплитуды, так как это 3D волна
            // Установка цвета в зависимости от амплитуды
            sf::Color color(static_cast<sf::Uint8>(amplitude * 255), 50, 50);
            image.setPixel(x, y, color);
        }
    }

    // Создаем текстуру изображения
    sf::Texture texture;
    texture.loadFromImage(image);

    // Создаем спрайт с текстурой
    sf::Sprite sprite(texture);

    // Отрисовываем спрайт
    window.draw(sprite);
}

int main() {
    // Создаем окно SFML
    sf::RenderWindow window(sf::VideoMode(800, 600), "3D Waves");

    // Пример трехмерного вектора амплитуд
    std::vector<std::vector<std::vector<float>>> amplitudes(800, std::vector<std::vector<float>>(600, std::vector<float>(1, 0.0f)));

    // Главный цикл программы
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Очистка экрана
        window.clear();

        // Отрисовка волн
        drawWaves(window, amplitudes);

        // Отображение отрисованного на экране
        window.display();
    }

    return 0;
}
