# Выбор компилятора
CC      = gcc

# Флаги компилятора: 
# -O2: Включить оптимизацию 
# -Wall: Показывать все предупреждения
# -I. (и другие): Где искать заголовочные файлы
CFLAGS  = -O2 -Wall -I. -Iinclude 

LDFLAGS = -lm

# Путь к заранее скомпилированной библиотеке CSPICE.
CSPICE  = lib/cspice.a

# Исходные файлы, которые нужно скомпилировать
SRC = dynamics.c integrator.c main.c
# Превращаем список .c файлов в список .o файлов
OBJ = $(SRC:.c=.o)

# Имя итогового исполняемого файла
TARGET = messenger_orbit


all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(CSPICE) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET) orbit_diff.txt