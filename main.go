package main

import (
    "bytes"
    "flag"
    "fmt"
    "errors"
    "strconv"
    "os/exec"
    "os"
    "math"
    "strings"
)

type Datetime struct {
    year int
    month int
    day int
    hour int
    minute int
}

type Location struct {
    latitude float64
    longitude float64
    standardMeridian int
}

func floatsToString(floats []float64) string {
    var builder strings.Builder
    for i, f := range floats {
        if i > 0 {
            builder.WriteRune(' ')
        }
        builder.WriteString(strconv.FormatFloat(f, 'f', -1, 64))
    }
    return builder.String()
}

func parseMonth(str string) (int, error) {
    month, err := strconv.Atoi(str)
    if err != nil {
        return 0, err
    }
    if (month < 1 || month > 12) {
        return 0, errors.New("Invalid month")
    }
    return month, nil
}

func parseDay(str string) (int, error) {
    day, err := strconv.Atoi(str)
    if err != nil {
        return 0, err
    }
    if (day < 1 || day > 31) {
        return 0, errors.New("Invalid day")
    }
    return day, nil
}

func parseHour(str string) (int, error) {
    hours, err := strconv.Atoi(str)
    if err != nil {
        return 0, err
    }
    if (hours < 0 || hours > 24) {
        return 0, errors.New("Invalid hours")
    }
    return hours, nil
}

func parseMinute(str string) (int, error) {
    minutes, err := strconv.Atoi(str)
    if err != nil {
        return 0, err
    }
    if (minutes < 0 || minutes >= 60) {
        return 0, errors.New("Invalid minutes")
    }
    return minutes, nil
}

func printVersion() {
    fmt.Println("genlrtsky v0.1")

    cmd := exec.Command("uvspec", "-v") 

    var stdoutBuf, stderrBuf bytes.Buffer
    cmd.Stdout = &stdoutBuf
    cmd.Stderr = &stderrBuf

    err := cmd.Run()
    if err != nil {
        panic(err)
    }
    fmt.Printf(stderrBuf.String())
    _, err = exec.LookPath("zenith")
    if err != nil {
        fmt.Println("zenith not found")
    }
    return
}

func generateUniformSamples(step int) ([]float64, []float64) {
    if 90 % step != 0 {
        panic("Angluar resolution not divisive by 90")
    }
    cos_thetas := make([]float64, 0, 90/step)
    for i := 0; i < 90; i += step {
        cos_thetas = append(cos_thetas, -math.Cos(float64(i) * (math.Pi / 180.0)))
    }
    phis := make([]float64, 0, 360/step)
    for i := 0; i < 360; i += step {
        phis = append(phis, float64(i))
    }
    return cos_thetas, phis
}

func getSolarAngle(dt Datetime, lc Location) (string, string) {
    cmd := exec.Command("zenith", "-y", strconv.Itoa(dt.year), 
        "-a", strconv.FormatFloat(lc.latitude, 'f', -1, 64),
        "-o", strconv.FormatFloat(lc.longitude, 'f', -1, 64),
        "-s", strconv.Itoa(lc.standardMeridian),
        strconv.Itoa(dt.day), strconv.Itoa(dt.month), 
        strconv.Itoa(dt.hour), strconv.Itoa(dt.minute),
    )
    var stdoutBuf, stderrBuf bytes.Buffer
    cmd.Stdout = &stdoutBuf
    cmd.Stderr = &stderrBuf

    fmt.Println(cmd.Args)
    err := cmd.Run()
    if err != nil {
        fmt.Fprintln(os.Stderr, stderrBuf.String())
        panic(err)
    }
    // stdout: 12:00:00   60.4563  358.4449
    solarAngles := strings.Fields(stdoutBuf.String())
    return solarAngles[1], solarAngles[2]
}


func main() {

    libPath := os.Getenv("LIBRADTRAN_DATA_FILES")
    if libPath == "" {
        fmt.Fprintln(os.Stderr, "LIBRADTRAN_DATA_FILES not set")
        return
    }

    // Command-line flags
    year := flag.Int("y", 2003, "Year")
    latitude := flag.Float64("a", 0.0, "latitude")
    longitude := flag.Float64("o", 0.0, "longitude")
    standardMeridian := flag.Int("m", 0, "Standard meridian, west positive")
    version := flag.Bool("version", false, "Print version")
    altitude := flag.Float64("altitude", 0.0, "Altitude")
    albedo := flag.Float64("albedo", 0.2, "Albedo")
    aerosol := flag.String("aerosol", "continental_average", "Standard aerosol profile name")
    step := flag.Int("step", 3, "Angular step")

    flag.Parse()

    if *version {
        printVersion()
    }

    // Check and get positional arguments before any optional flags.
    positionalArgs := []string{}
    for _, arg := range os.Args[1:] { // Skipping the program name.
        if strings.HasPrefix(arg, "-") {
            break // Stop on the first flag.
        }
        positionalArgs = append(positionalArgs, arg)
    }

    // Make sure we got the first three positional arguments.
    if len(positionalArgs) < 4 {
        fmt.Fprintf(os.Stderr, "Usage: %s month day hour minute [options]\n", os.Args[0])
        os.Exit(1)
    }

    month, err := parseMonth(positionalArgs[0])
    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        return
    }
    day, err := parseDay(positionalArgs[1])
    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        return
    }
    hour, err := parseHour(positionalArgs[2])
    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        return
    }
    minute, err := parseMinute(positionalArgs[3])
    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        return
    }


    dt := Datetime{year: *year, month: month, day: day, hour: hour, minute: minute}
    lc := Location{latitude: *latitude, longitude: *longitude, standardMeridian: *standardMeridian}
    // Get the solar zenith angle
    solarZenith, solarAzimuth := getSolarAngle(dt, lc)

    umu, phis := generateUniformSamples(*step)

    input := "data_files_path " + libPath + "\n"
    input += "source solar apm_1nm\n"
    input += "pseudospherical\n"   // use pseudospherical solver
    input += "altitude " + strconv.FormatFloat(*altitude, 'f', -1, 64) + "\n"
    input += "albedo " + strconv.FormatFloat(*albedo, 'f', -1, 64) + "\n"
    input += "sza " + solarZenith + "\n"
    input += "phi0 " + solarAzimuth + "\n"
    input += "aerosol_default\n"
    input += "aerosol_species_file " + *aerosol + "\n"
    input += "umu " + floatsToString(umu) + "\n"
    input += "phi " + floatsToString(phis) + "\n"
    input += "wavelength 550\n"
    input += "output_user lambda edir edn uu\n"

    fmt.Println(input)
    cmd := exec.Command("uvspec")
    cmd.Stdin = bytes.NewBufferString(input)

    // Capture the output
    var stderr, out bytes.Buffer
    cmd.Stdout = &out
    cmd.Stderr = &stderr

    fmt.Println(cmd.Args)
    // Start and wait for the command to finish
    if err := cmd.Run(); err != nil {
        fmt.Println(stderr.String())
        panic(err)
    }

    // Print the output from the command
    fmt.Println("Command output:", out.String())

}
