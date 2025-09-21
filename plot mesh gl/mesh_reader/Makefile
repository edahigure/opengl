CC = gcc
CFLAGS = -Wall -O2
TARGET = mesh_reader
SOURCES = main.c mesh_reader.c

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCES)

clean:
	rm -f $(TARGET) *.o

.PHONY: clean
