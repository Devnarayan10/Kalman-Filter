import serial
import csv
import time

# Set the serial port and baud rate
serial_port = 'COM4'  # Change this to your Arduino's COM port
baud_rate = 9600
file_name = 'data.csv'

# Create a serial connection
ser = serial.Serial(serial_port, baud_rate)
time.sleep(1)  # Wait for the connection to be established

print("Serial connection established.")

# Open a CSV file for writing
with open(file_name, mode='a', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)

    try:
        while True:
            # Read a line from the serial port
            try:
                line = ser.readline().decode('utf-8').rstrip()
                if line:
                    print("Received data:", line)
                    # Split the line into a list
                    data = line.split(',')
                    # Write the data to the CSV file
                    csv_writer.writerow(data)
                    csv_file.flush()
                    print("Data written to CSV file.")
                else:
                    print("No data received.")
            except serial.SerialException as e:
                print("Serial error:", e)
    except KeyboardInterrupt:
        print("Logging stopped.")
    finally:
        ser.close()
        print("Serial connection closed.")