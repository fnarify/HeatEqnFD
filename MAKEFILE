CC=gcc
CFLAGS=-Wall -Wextra -g
TARGET=heat_eqn

.PHONY: heat_eqn clean

heat_eqn:
	$(CC) -o $(TARGET) $(CFLAGS) $(TARGET).c

clean:
	rm -f $(TARGET)
