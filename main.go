package main

import (
	"bytes"
	"encoding/binary"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"text/template"
	"sync"
)

// Manuanlly defined constants
const NSSAMP = 20
const WVLSPAN = 400  // Wavelength span in nm 780 - 380
const XRES = 32
const YRES = 32

type InputParams struct {
	LibPath            string
	SolFile            string
	GroundAlbedo       float64
	SolarZenith        string
	SolarAzimuth       string
	AerosolProfileName string
	CosTheta           float64
	Phi                float64
	WavelengthGridFile string
	Output             string
}

type Datetime struct {
	year   int
	month  int
	day    int
	hour   int
	minute int
}

type Location struct {
	latitude         float64
	longitude        float64
	standardMeridian int
}

type Result struct {
        X, Y int
        Data [NSSAMP+1]uint8 // Data to be written to hsrf
	Err  error
}


func panicError(e error) {
	if e != nil {
		panic(e)
	}
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
	if month < 1 || month > 12 {
		return 0, errors.New("Invalid month")
	}
	return month, nil
}

func parseDay(str string) (int, error) {
	day, err := strconv.Atoi(str)
	if err != nil {
		return 0, err
	}
	if day < 1 || day > 31 {
		return 0, errors.New("Invalid day")
	}
	return day, nil
}

func parseHour(str string) (int, error) {
	hours, err := strconv.Atoi(str)
	if err != nil {
		return 0, err
	}
	if hours < 0 || hours > 24 {
		return 0, errors.New("Invalid hours")
	}
	return hours, nil
}

func parseMinute(str string) (int, error) {
	minutes, err := strconv.Atoi(str)
	if err != nil {
		return 0, err
	}
	if minutes < 0 || minutes >= 60 {
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
	panicError(err)
	fmt.Printf(stderrBuf.String())
	_, err = exec.LookPath("zenith")
	if err != nil {
		fmt.Println("zenith not found")
	}
	return
}

func generateUniformSamples(step int) ([]float64, []float64) {
	if 90%step != 0 {
		panic("Angluar resolution not divisive by 90")
	}
	cos_thetas := make([]float64, 0, 90/step)
	for i := 0; i < 90; i += step {
		cos_thetas = append(cos_thetas, -math.Cos(float64(i)*(math.Pi/180.0)))
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

	// fmt.Println(cmd.Args)
	err := cmd.Run()
	panicError(err)
	// stdout: 12:00:00   60.4563  358.4449
	solarAngles := strings.Fields(stdoutBuf.String())
	return solarAngles[1], solarAngles[2]
}

func pix2loc(xres int, yres int, px int, py int) (float64, float64) {
	// Compute image location from pixel posistion
	// Assuming YMAJOR and YDECR
	x := px
	y := py
	y = yres - 1 - y
	loc0 := (float64(x) + 0.5) / float64(xres)
	loc1 := (float64(y) + 0.5) / float64(yres)
	return loc0, loc1
}

func viewray(x float64, y float64) (float64, float64) {
	vdir := [3]float64{0.0, 0.0, 1.0}
	hvec := [3]float64{-1.0, 0.0, 0.0}
	vvec := [3]float64{0.0, 1.0, 0.0}
	x -= 0.5
	y -= 0.5
	d := math.Sqrt(x*x + y*y)
	z := math.Cos(math.Pi * d)
	if d < 1e-6 {
		d = math.Pi
	} else {
		d = math.Sqrt(1-z*z) / d
	}
	x *= d
	y *= d
	direc := [3]float64{
		z*vdir[0] + x*hvec[0] + y*vvec[0],
		z*vdir[1] + x*hvec[1] + y*vvec[1],
		z*vdir[2] + x*hvec[2] + y*vvec[2],
	}
	phi := math.Atan2(direc[1], direc[0]) * 180 / math.Pi
	return -direc[2], phi
}

func scolor2scolr(scol [20]float64, ncs int) [21]uint8 {
	const COLXS uint8 = 128
        var p float64
	sclr := [21]uint8{}

        // Find the largest value in scol
        p = scol[0]
        for i := 1; i < ncs; i++ {
                if scol[i] > p {
                        p = scol[i]
                }
        }
        if p > 1e-32 {
                frac, exp := math.Frexp(p)
                p = frac * 256.0 / p // Ensure this isn't zero by checking `p` and `frac`
                sclr[ncs] = uint8(exp) + COLXS
                for i := 0; i < ncs; i++ {
                        if scol[i] > 0 {
                                sclr[i] = uint8(scol[i] * p)
                        }
                }
        } else {
                // Fill sclr with zeros if max value is not significant
                for i := range sclr {
                        sclr[i] = 0
                }
        }
	return sclr
}

// putbinary is a function to write small objects to an io.Writer like a file.
func putbinary(w io.Writer, data interface{}) (int, error) {
        buf := new(bytes.Buffer)
        if err := binary.Write(buf, binary.LittleEndian, data); err != nil {
                return 0, err
        }
        bytes := buf.Bytes()

        // Use io.Writer's Write method directly if large object
        if len(bytes) > 128 {
                return w.Write(bytes)
        }

        // Write byte-by-byte for small objects
        count := 0
        for _, b := range bytes {
                if _, err := w.Write([]byte{b}); err != nil {
                        return count, err
                }
                count++
        }
        return count, nil
}

func worker(x, y int, input InputParams, tmpl *template.Template, wg *sync.WaitGroup, concurrencyLimits chan struct{}, results chan<- Result) {
        // Perform CPU-heavy `uvspec` work...
        // (Your operations here)
	defer wg.Done()
	defer func() { <-concurrencyLimits }() // Release the "slot"

        // Simulated result and error
        var result Result

        result.X = x
        result.Y = y

        // Imagine `myProcessingFunc` is your processing function that calls `uvspec`
        result.Data, result.Err = compute(x, y, input, tmpl)

        // Send the result to the results channel
        results <- result
}

func compute(x, y int, input InputParams, tmpl *template.Template) ([NSSAMP+1]uint8, error) {

	loc0, loc1 := pix2loc(XRES, YRES, x, y)
	umu, phi := viewray(loc0, loc1)

	var templateBuffer bytes.Buffer
	input.CosTheta = umu
	input.Phi = phi
	err := tmpl.Execute(&templateBuffer, input)
	if err != nil {
		panic(err)
	}
	inputStr := templateBuffer.String()

	tdir, err := os.MkdirTemp("", "")
	if err != nil {
		panic(err)
	}
	defer os.RemoveAll(tdir)

	// fmt.Println(inputStr)
	cmd := exec.Command("uvspec")
	cmd.Stdin = bytes.NewBufferString(inputStr)
	cmd.Dir = tdir

	// Capture the output
	var stderr, out bytes.Buffer
	cmd.Stderr = &stderr
	cmd.Stdout = &out

	// fmt.Println(cmd.Args)
	err = cmd.Run()
	// fmt.Fprintln(os.Stderr, stderr.String())
	if err != nil {
		panic(err)
	}
	res, err := os.ReadFile(filepath.Join(tdir, "mc.rad.spc"))
	if err != nil {
		panic(err)
	}
	rows := strings.Split(string(res), "\n")
	scolor := [20]float64{}
	for i, row := range rows {
		if row == "" {
			continue
		}
		fields := strings.Fields(row)
		value, err := strconv.ParseFloat(fields[4], 64)
		if err != nil {
			panic(err)
		}
		scolor[i] = value * WVLSPAN / 1000  // mW to W
	}
	sclr := scolor2scolr(scolor, NSSAMP)
	fmt.Println(x, y)
	return sclr, err
}

func main() {


	inputTemplate := `data_files_path {{.LibPath}} 
source solar {{.SolFile}}
rte_solver mystic
albedo {{.GroundAlbedo}}
sza {{.SolarZenith}}
phi0 {{.SolarAzimuth}}
aerosol_default
aerosol_species_file {{.AerosolProfileName}}
umu {{.CosTheta}}
phi {{.Phi}}
wavelength_grid_file {{.WavelengthGridFile}}
mol_abs_param crs
mc_vroom on`

	tmpl, err := template.New("input").Parse(inputTemplate)
	if err != nil {
		panic(err)
	}

	input := InputParams{}

	// Extraterrestrial solar irradiance in mW/m^2/nm
	EXTSOL := [NSSAMP]float64{1192.74, 1772.2, 1394.52, 2069.67,
		2074.3, 1964.5, 1971.29, 1942.59, 1870.53,
		1852.63, 1772.84, 1737.29, 1662.2, 1586.51,
		1537.19, 1475.65, 1400.21, 1320.18, 1273.23,
		1220.45,
	}

	libPath := os.Getenv("LIBRADTRAN_DATA_FILES")
	if libPath == "" {
		fmt.Fprintln(os.Stderr, "LIBRADTRAN_DATA_FILES not set")
		return
	}
	input.LibPath = libPath

	// Command-line flags
	year := flag.Int("y", 2003, "Year")
	latitude := flag.Float64("a", 37.0, "latitude")
	longitude := flag.Float64("o", 122.0, "longitude")
	standardMeridian := flag.Int("m", 120, "Standard meridian, west positive")
	version := flag.Bool("version", false, "Print version")
	// altitude := flag.Float64("altitude", 0.0, "Altitude")
	albedo := flag.Float64("albedo", 0.2, "Albedo")
	aerosol := flag.String("aerosol", "continental_average", "Standard aerosol profile name")
	// cloudCover := flag.String("cloud cover", "continental_average", "Standard aerosol profile name")
	// step := flag.Int("step", 3, "Angular step")

	input.GroundAlbedo = *albedo
	input.AerosolProfileName = *aerosol

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
		fmt.Fprintf(os.Stderr, "Usage: genlrtsky month day hour minute [options]\n")
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

	// Get the solar angle
	dt := Datetime{year: *year, month: month, day: day, hour: hour, minute: minute}
	lc := Location{latitude: *latitude, longitude: *longitude, standardMeridian: *standardMeridian}
	solarZenith, solarAzimuth := getSolarAngle(dt, lc)
	input.SolarZenith = solarZenith
	input.SolarAzimuth = solarAzimuth

	// umu, phis := generateUniformSamples(*step)

	// Create a temp file to store solar and wavelength data
	wavelengthGridFile, err := os.CreateTemp("", "lambda.dat")
	panicError(err)
	defer os.Remove(wavelengthGridFile.Name())
	defer wavelengthGridFile.Close()
	input.WavelengthGridFile = wavelengthGridFile.Name()

	solFile, err := os.CreateTemp("", "solar.dat")
	if err != nil {
		panic(err)
	}
	defer os.Remove(solFile.Name())
	defer solFile.Close()
	input.SolFile = solFile.Name()

	for index, value := range EXTSOL {
		wavelength := 390 + index*NSSAMP
		wavelengthStr := strconv.Itoa(wavelength)
		wavelengthGridFile.WriteString(wavelengthStr + "\n")
		solFile.WriteString(wavelengthStr + " " + strconv.FormatFloat(value, 'f', -1, 64) + "\n")
	}
	results := make(chan Result, YRES*XRES)
	wg := &sync.WaitGroup{}

	// Limit the number of concurrent goroutines to the number of processor cores
        numCores := runtime.NumCPU()
	// numCores = 1
	// numCores *= 2
        concurrencyLimit := make(chan struct{}, numCores) // Semaphore-like channel

        for y := 0; y < YRES; y++ {
                for x := 0; x < XRES; x++ {
                        concurrencyLimit <- struct{}{} // Acquire a "slot"
                        wg.Add(1)
                        go worker(x, y, input, tmpl, wg, concurrencyLimit, results)
                }
        }

        // Ensure all workers are finished before closing the results channel
        go func() {
                wg.Wait()
                close(results)
        }()

	fmt.Println("Done with workers...")

        orderedResults := make([][]Result, YRES)
        for y := range orderedResults {
                orderedResults[y] = make([]Result, XRES)
        }
	
	fmt.Println("Assigning resutls...")

        for res := range results {
		fmt.Println(res.X, res.Y)
                if res.Err != nil {
                        fmt.Printf("Error processing (%d, %d): %v\n", res.X, res.Y, res.Err)
                        continue
                }
                orderedResults[res.Y][res.X] = res
        }

	fmt.Println("Writing results...")

        // Assuming you want to accumulate all outputs into `hsrf`
        var hsr bytes.Buffer

        // Write the results in `hsrf` in order
        for y := range orderedResults {
                for x := range orderedResults[y] {
			fmt.Println(x, y)
                        if res := orderedResults[y][x]; res.Err == nil {
                                for _, value := range res.Data {
                                        if err := binary.Write(&hsr, binary.LittleEndian, value); err != nil {
                                                panic(err)
                                        }
                                }
                        }
                }
        }

        // Now `hsrf` contains the ordered results
        // (You'll want to write these to a file or process further as needed)

	// hsrf, err := os.OpenFile("sky.hsr", os.O_CREATE|os.O_WRONLY, 0644)
	hsrf, err := os.Create("sky.hsr")
	defer hsrf.Close()
	hsrf.WriteString("#?RADIANCE\n")
	hsrf.WriteString("NCOMP=20\n")
	hsrf.WriteString("WAVELENGTH_SPLITS= 380 480 588 780\n")
	hsrf.WriteString("FORMAT=Radiance_spectra\n\n")
	resStr := fmt.Sprintf("-Y       %d +X       %d\n", YRES, XRES)
	hsrf.WriteString(resStr)
	hsrf.Write(hsr.Bytes())

	// for y := 0; y < YRES; y++ {
	// 	for x := 0; x < XRES; x++ {
	// 		fmt.Println(x, y)
	// 		loc0, loc1 := pix2loc(XRES, YRES, x, y)
	// 		umu, phi := viewray(loc0, loc1)
	//
	// 		var templateBuffer bytes.Buffer
	// 		input.CosTheta = umu
	// 		input.Phi = phi
	// 		err = tmpl.Execute(&templateBuffer, input)
	// 		if err != nil {
	// 			panic(err)
	// 		}
	// 		inputStr := templateBuffer.String()
	//
	// 		tdir, err := os.MkdirTemp("", "")
	// 		if err != nil {
	// 			panic(err)
	// 		}
	// 		defer os.RemoveAll(tdir)
	//
	// 		// fmt.Println(inputStr)
	// 		cmd := exec.Command("uvspec")
	// 		cmd.Stdin = bytes.NewBufferString(inputStr)
	// 		cmd.Dir = tdir
	//
	// 		// Capture the output
	// 		var stderr, out bytes.Buffer
	// 		cmd.Stderr = &stderr
	// 		cmd.Stdout = &out
	//
	// 		// fmt.Println(cmd.Args)
	// 		err = cmd.Run()
	// 		// fmt.Fprintln(os.Stderr, stderr.String())
	// 		if err != nil {
	// 			fmt.Fprintln(os.Stderr, stderr.String())
	// 			return
	// 		}
	// 		res, err := os.ReadFile(filepath.Join(tdir, "mc.rad.spc"))
	// 		if err != nil {
	// 			panic(err)
	// 		}
	// 		rows := strings.Split(string(res), "\n")
	// 		scolor := [20]float64{}
	// 		for i, row := range rows {
	// 			if row == "" {
	// 				continue
	// 			}
	// 			fields := strings.Fields(row)
	// 			value, err := strconv.ParseFloat(fields[4], 64)
	// 			if err != nil {
	// 				panic(err)
	// 			}
	// 			scolor[i] = value * WVLSPAN / 1000  // mW to W
	// 		}
	// 		sclr := scolor2scolr(scolor, NSSAMP)
	// 		for _, value := range sclr {
	// 			err = binary.Write(hsrf, binary.LittleEndian, value)
	// 			if err != nil {
	// 				panic(err)
	// 			}
	// 		}
	// 	}
	// }

	// dni := [NSSAMP]string{}
	// uu := [NSSAMP * 121 * 31]string{}
	// rows := strings.Split(out.String(), "\n")
	// for index, row := range rows {
	// 	if row == "" {
	// 		continue
	// 	}
	// 	fields := strings.Fields(row)
	// 	fmt.Println(len(fields))
	// 	dni[index] = fields[1]
	// 	for index2, value2 := range fields[3:] {
	// 		idx := index2 * NSSAMP + index
	// 		uu[idx] = value2
	// 	}
	// }
	//
	// datFile, err := os.Create("res.dat")
	// if err != nil {
	// 	panic(err)
	// }
	// defer datFile.Close()
	// datFile.WriteString("3\n")
	// datFile.WriteString("0 360 121\n")
	// datFile.WriteString("0 90 31\n")
	// datFile.WriteString("380 780 20\n")
	// for _, value := range uu {
	// 	datFile.WriteString(value + "\n")
	// }

	// numWavelengths := 20
	// resultsChan := make(chan string, numWavelengths) // Buffered channel
	// var wg sync.WaitGroup

	// for i := 390; i <= 770; i += 20 {
	// 	// wg.Add(1)
	//
	// 	// go func(wavelength int) {
	// 		defer wg.Done()
	// 		input += "wavelength " + strconv.Itoa(wavelength) + "\n"
	// 		fmt.Println(input)
	// 		cmd := exec.Command("uvspec")
	// 		cmd.Stdin = bytes.NewBufferString(input)
	//
	// 		// Capture the output
	// 		var stderr, out bytes.Buffer
	// 		cmd.Stdout = &out
	// 		cmd.Stderr = &stderr
	//
	// 		fmt.Println(cmd.Args)
	// 		// Start and wait for the command to finish
	// 		err := cmd.Run()
	// 		if err != nil {
	// 			fmt.Fprintln(os.Stderr, stderr.String())
	// 			panic(err)
	// 		}
	// 		resultsChan <- out.String()
	// 	}(i)
	//
	// 	// Print the output from the command
	// 	// fmt.Println("Command output:", out.String())
	// }

	// Start a collector goroutine
	// go func() {
	// 	wg.Wait()
	// 	close(resultsChan)
	// }()

	// Collect the ordered results
	// var results []string
	// for r := range resultsChan {
	// 	nums := strings.Fields(r)
	// 	fmt.Println(nums[:8])
	// 	// results = append(results, r)
	// }

	// fmt.Println("Results:")
	// for _, r := range results {
	// 	fmt.Println(r)
	// }

}
